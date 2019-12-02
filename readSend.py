import socket
import array
import time

import numpy as np
import sys

# Create a TCP/IP socket
sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

# Bind the socket to the port
portAddressServer = 12000
portAddressClient = 12001
clientAddress = '10.0.10.100'
serverAddress = '10.0.10.101'
numberSamples = 2000
numberValueFromTank = 6
packetSize = 48
server_address = (serverAddress, portAddressServer)
print('starting up on {} port {}\n\n'.format(server_address, portAddressServer))
sock.bind(server_address)

dataSystem = np.zeros(shape=[numberSamples,numberValueFromTank],dtype=float)
index = 0
running = True

while running:

    #print('\nwaiting to receive message')
    data, address = sock.recvfrom(packetSize)
    #print('\nreceived {} bytes from {}'.format(len(data), serverAddress))
    doubles_sequence = array.array('d', data)
    #doubles_sequence.byteswap() #default is little indian, if big indian then uncomment this
    print('{},{},{},{},{},{}'.format(doubles_sequence[0], doubles_sequence[1], doubles_sequence[2],doubles_sequence[3],doubles_sequence[4],doubles_sequence[5]))
    dataSystem[index] = doubles_sequence[0], doubles_sequence[1], doubles_sequence[2],doubles_sequence[3],doubles_sequence[4],doubles_sequence[5]
    index +=1

    if index >= numberSamples:
        running = False;
    #if data:
    #   sent = sock.sendto(data, address)
    #//   print('\nsent {} bytes back to {}'.format(len(sent), clientAddress))