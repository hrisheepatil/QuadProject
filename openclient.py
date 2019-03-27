import io
import socket
import struct
from PIL import Image
import cv2
import numpy as np

# Start a socket listening for connections on 0.0.0.0:8000 (0.0.0.0 means
# all interfaces)
server_socket = socket.socket()
server_socket.bind(('192.168.1.3', 8000))
server_socket.listen(0)

# Accept a single connection and make a file-like object out of it
connection = server_socket.accept()[0].makefile('rb')
try:
    while True:
        # Read the length of the image as a 32-bit unsigned int. If the
        # length is zero, quit the loop
        image_len = struct.unpack('<L', connection.read(struct.calcsize('<L')))[0]
        if not image_len:
            break
        image_stream = io.BytesIO()
        image_stream.write(connection.read(image_len))
        image_stream.seek(0)
        image = Image.open(image_stream)
        print(image)
        cv_image = np.array(image)
        print(cv_image)
        cv2.imshow('cv_image',cv_image)
        if cv2.waitKey(1) & 0xFF == ord('q'):
            break

finally:
    connection.close()
    server_socket.close()
