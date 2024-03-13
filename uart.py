#pip install pyserial
import serial

file_path = '/Users/luana/Desktop/tiempos.txt'  # Ruta del archivo donde se guardarán los datos

# Configura la conexión serial con los parámetros apropiados
ser = serial.Serial('/dev/cu.usbserial-0001', 57600)

try:
    with open(file_path, 'a') as file:
        #print("dt, x, y, z, vx, vy, vz, q1, q2, q3, q4, traza_cov")
        #file.write("dt, x, y, z, vx, vy, vz, q1, q2, q3, q4, traza_cov\n")
        while True:
            # Lee datos del dispositivo UART
            data = ser.readline().decode().strip()
            print(data)
            # Filtra líneas vacías o que contienen solo espacios en blanco
            if data and not data.isspace():
                file.write(data + '\n') # Escribe los datos en el archivo
except KeyboardInterrupt:
    ser.close()  # Cierra la conexión serial al interrumpir el programa con Ctrl+C

