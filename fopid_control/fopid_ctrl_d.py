from fopid_control import fpid
from struct import pack, unpack
import socket
from select import select
import time

'''
This is a proof-of-concept code for running real-time
FO PID control from Python using a UDP socket
'''

UDP_IP = "127.0.0.1"
UDP_PORT_REMOTE = 5005  # The port where the computed control law is sent (client port)
UDP_PORT_LOCAL = 5006  # Port to which the error signal must be sent (server port)
UDP_PORT_CTRL = 5201  # Control port. Allows to, e.g., set controller parameters (server port)
UDP_BUFFER_SIZE = 4096

# Set the initial FOPID parameters. The parameters of the approximation are NOT set in this example
params = {"fpid": {"Kp": 1, "Ki": 1, "Kd": 1, "lam": 0.5, "mu": 0.5}}

# This initializes and computes the controller. From now on, you can access all of the data
# of this approximation *and* can run the control algorithm.
fpid_c = fpid.fpid_2iir(params)

# Latest computed control law
last_out = 0

# Change the parameters of the FOPID controller
def change_serv_params(ctrlsock):
    global fpid_c
    data, addr = ctrlsock.recvfrom(UDP_BUFFER_SIZE)
    inp = unpack("<cddddd", data)  # < = little endian, c = character (command), d = doubles

    # What is the command?
    command = inp[0].decode('ascii')

    # c: change FOPID parameters
    if command == "c":
        params = {"fpid": {"Kp": inp[1], "Ki": inp[2], "Kd": inp[3], "lam": inp[4], "mu": inp[5]}}
        fpid_c = fpid.fpid_2iir(params)
        print("FOPID parameter change requested:")
        print(params)
        print("New controller parameters were successfully applied.")


# Compute the control law
def do_control_io(locsock, remsock):

    # Global function in lieu of class definition
    global last_out, fpid_c
    currentTimeIn = time.time()
    # Receive the data from the socket
    data, addr = locsock.recvfrom(UDP_BUFFER_SIZE)

    # Decode data
    inp = unpack("<d", data)[0]  # < = little endian, d = double

    # Run the control algorithm
    if not fpid_c.COMPUTING_FOPID_APPROXIMATION:
        out = fpid_c.compute_fopid_control_law(inp)
        last_out = out
    else:
        out = last_out

    # Send the control law
    msg = pack("<d", out)
    remsock.sendto(msg, (UDP_IP, UDP_PORT_REMOTE))
    currentTimeOut = time.time()
    print("Recieved:{0:.4f},\t\t sent:{2:.4f}, \t\t dtCompute = {4:.5f}".format(inp, data, out, msg,
                                                                                currentTimeOut - currentTimeIn))

def startControl():
    print("Starting control system server...")
    print("Server IP:", UDP_IP)
    print("Server port [receive e(t) from controlled process]:", UDP_PORT_LOCAL)
    print("Server control port [receive modified FOPID parameters]:", UDP_PORT_CTRL)
    print("Remote (client) port [send computed control law u(t)]:", UDP_PORT_REMOTE)
    # Local socket: server
    locsock = socket.socket(socket.AF_INET,  # Internet
                            socket.SOCK_DGRAM)  # UDP
    # Remote socket: client
    remsock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    # Local socket for controlling the server: server
    ctrlsock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    locsock.bind((UDP_IP, UDP_PORT_LOCAL))
    ctrlsock.bind((UDP_IP, UDP_PORT_CTRL))
    inp = [locsock, ctrlsock]
    print("Server started successfully. Waiting for communications.")
    # Run the select loop
    while True:

        iready, _, _ = select(inp, [], [])
        for s in iready:
            if s == locsock:
                do_control_io(locsock, remsock)
            elif s == ctrlsock:
                change_serv_params(ctrlsock)

if __name__ == '__main__':
    startControl()