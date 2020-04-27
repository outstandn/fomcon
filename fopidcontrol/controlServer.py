from fopidcontrol import fpid
from struct import pack, unpack
import socket
from select import select
import time
import sys
from multiprocessing import Pipe
import traceback
import ast
import asyncio
from addict import Dict


DURATION = 600
SLEEPTIME = 20
shareMemoGUI, shareMemoControl = Pipe()
SERVERRUNNIN =False

class fomconControlServer():
    '''
    This is a proof-of-concept code for running real-time
    FO PID control from Python using a UDP socket
    '''

    def __init__(self,time2Run,lock =None):
        self.UDP_IP_RECV_LOCAL = "127.0.0.1"
        self.UDP_IP_SEND2REMOTE = "127.0.0.1"
        self.UDP_PORT_CTRL = 5201               # Port for control Should be set manaually in here..and should be same in GUI
        self.UDP_PORT_REMOTE_IN = 5005          # The port where the computed control law is sent (client port)
        self.UDP_PORT_LOCAL = 5006              # Port to which the error signal must be sent (server port)
        self.UDP_BUFFER_SIZE = 6144             # Control port. Allows to, e.g., set controller parameters (server port)
        self.time2Run = time2Run
        self.lock = lock
        self.start = False
        self.exitt = False
        self.plantIPPortConnected = False

        # Local socket: server
        self.locsock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)  # Internet  # UDP

        # Remote socket: client
        self.remsock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

        # Local socket for controlling the server: server
        self.ctrlsock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        # Bind Receiveing
        # self.locsock.bind((self.UDP_IP_RECV_LOCAL, self.UDP_PORT_LOCAL))
        self.ctrlsock.bind((self.UDP_IP_RECV_LOCAL, self.UDP_PORT_CTRL))

        # Set the initial FOPID parameters. The parameters of the approximation are NOT set in this example
        self.params = {"fpid": {"Kp": 1, "Ki": 1, "Kd": 1, "lam": 0.5, "mu": 0.5}}

        # This initializes and computes the controller. From now on, you can access all of the data
        # of this approximation *and* can run the control algorithm.
        self.fpid_c = fpid.fpid_2iir(self.params)

        # Latest computed control law
        self.last_out = 0

    def updateBindingOnly(self):
        try:
            #Bind Receiveing
            self.locsock.bind((self.UDP_IP_RECV_LOCAL, self.UDP_PORT_LOCAL))
            return True
        except:
            return False

    def disConnectAndUpdate(self, ipRecv, recvPort, ipSend, sendPort,simetime):
        try:
            # disconnect
            self.closePlantCon()
            time.sleep(SLEEPTIME)
            #update IP and Port
            self.UDP_IP_RECV_LOCAL, self.UDP_PORT_LOCAL = ipRecv, recvPort
            self.UDP_IP_SEND2REMOTE,self.UDP_IP_SEND2REMOTE = ipSend, sendPort
            self.time2Run = simetime

            #Bind Receiveing
            self.locsock.bind((self.UDP_IP_RECV_LOCAL, self.UDP_PORT_LOCAL))

            return True
        except Exception:
            traceback.print_exc()
            return False

    #TODO: Use shared memeory to update IP and PORT or VIA UDP SOCKET
    def UpdateSocksMuilti(self,sharedMem):
        #disconnect
        self.closePlantCon()
        #update
        self.UDP_IP_RECV_LOCAL = sharedMem[0]
        self.UDP_IP_SEND2REMOTE = sharedMem[1]
        self.UDP_PORT_REMOTE_IN = sharedMem[3]
        self.UDP_PORT_LOCAL = sharedMem[2]
        self.UDP_PORT_CTRL = sharedMem[4]

    def closePlantCon(self):
        self.locsock.shutdown(socket.SHUT_RDWR)
        # self.remsock.shutdown(socket.SHUT_RDWR)
        self.locsock.close()
        # self.remsock.close()

        # Change the parameters of the FOPID controller
    def change_serv_params(self):
        data, addr = self.ctrlsock.recvfrom(self.UDP_BUFFER_SIZE)
        inp = data.decode('utf-8')
        input = ast.literal_eval(inp)
        # What is the command?
        command = list(input.keys())[0]

        if command == 'control':
            inputDict = Dict(input['control'])
            self.params = dict(fpid = dict(Kp=float(inputDict.Kp), Ki = float(inputDict.Ki), Kd = float(inputDict.Kd),lam = float(inputDict.lam), mu = float(inputDict.mu)))
            self.fpid_c = fpid.fpid_2iir(self.params)
            print("FOPID parameter change requested:")
            print(self.params)
            print("New controller parameters were successfully applied.")
        elif command == 'time':
            inputDict = Dict(input['time'])
            # time: Update simulation time
            self.time2Run = int(inputDict.time2Run)

        elif command == 'start':
            # start: Start controller FOPID parameters
            try:
                inputDict = Dict(input['start'])
                if self.plantIPPortConnected:
                    self.disConnectAndUpdate(inputDict.recvIP,int(inputDict.recvPORT),inputDict.sendIP,int(inputDict.sendPORT),int(inputDict.time2Run))
                else:
                    self.UDP_IP_SEND2REMOTE = inputDict.sendIP
                    self.UDP_IP_RECV_LOCAL = inputDict.recvIP
                    self.UDP_PORT_REMOTE_IN = int(inputDict.sendPORT)
                    self.UDP_PORT_LOCAL = int(inputDict.recvPORT)
                    self.time2Run = int(inputDict.time2Run)
                    self.updateBindingOnly()
                    self.start = True
                    print("Servers' Control port [receive modified FOPID parameters]:", self.UDP_PORT_CTRL)
                    print("Servers' Receive IP:Port = ", self.UDP_IP_RECV_LOCAL, ':', self.UDP_PORT_LOCAL)
                    print("Remote (client) IP:Port= ", self.UDP_IP_SEND2REMOTE, ':', self.UDP_PORT_REMOTE_IN)
                self.plantIPPortConnected = True
            except:
                self.plantIPPortConnected = False
                traceback.print_exc()
        elif command == 'stop':
            # stop: Stop controller FOPID parameters
            self.start = False
            self.time2Run = 0
            self.last_out = 0

        elif command == 'exit':
            # exit: Stop controller FOPID parameters
            self.exitt = True
            self.start = False
            self.time2Run = 0

    # Compute the control law
    def do_control_io(self):

        # Global function in lieu of class definition
        currentTimeIn = time.time()
        # Receive the data from the socket
        data, addr = self.locsock.recvfrom(self.UDP_PORT_CTRL)

        # Decode data
        inp = unpack("<d", data)[0]  # < = little endian, d = double

        # Run the control algorithm
        if not self.fpid_c.COMPUTING_FOPID_APPROXIMATION:
            out = self.fpid_c.compute_fopid_control_law(inp)
            self.last_out = out
        else:
            out = self.last_out

        # Send the control law
        msg = pack("<d", out)
        self.remsock.sendto(msg, (self.UDP_IP_SEND2REMOTE, self.UDP_PORT_REMOTE_IN))
        currentTimeOut = time.time()
        print("Recieved:{0:.4f},\t\t sent:{2:.4f}, \t\t dtCompute = {4:.5f}".format(inp, data, out, msg,
                                                                                    currentTimeOut - currentTimeIn))

    def startControl(self,timing=None):
        print("Waiting For Plant IP and Ports...")
        # self.updateBindingOnly()
        inp = [self.locsock,self.ctrlsock]
        # print("Server started successfully. Waiting for communications.")

        # Run the select loop
        self.time2Run = DURATION if timing is None else timing
        while self.exitt == False:
            t0 = time.time_ns()/1e9
            iready, _, _ = select(inp, [], [])
            for s in iready:
                if s == self.ctrlsock:
                    self.change_serv_params()

            while self.start:
                iready, _, _ = select(inp, [], [])
                for s in iready:
                    if s == self.locsock:
                        self.do_control_io()
                    elif s == self.ctrlsock:
                        self.change_serv_params()

                if ((time.time_ns()/1e9) - t0) >= self.time2Run:
                    self.start = False
                    break

        self.closePlantCon()
        self.ctrlsock.shutdown(socket.SHUT_RDWR)
        self.ctrlsock.close()

if __name__ == '__main__':
    timeargv = sys.argv[1] if len(sys.argv) == 2 else DURATION
    conn = fomconControlServer(timeargv, None)
    conn.startControl()

    #function for if called by a process. Time to run is in seconds (int),
    # sharedMem is a multiprocessing.Quene,use Queue.put and Queue.get
def controlServerAutoStart(dur, lock):
    con = fomconControlServer(dur, lock)
    con.startControl()
    time.sleep(5)