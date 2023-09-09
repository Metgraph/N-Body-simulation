import cv2

# Impostazioni del video
output_file = 'output_video.avi'
frame_rate = 150
frame_size = (1024//2, 768//2)

image_files = [f'images/2022-{i}.png' for i in range(1000)]

# Crea l'oggetto VideoWriter
fourcc = cv2.VideoWriter_fourcc(*'XVID')
out = cv2.VideoWriter(output_file, fourcc, frame_rate, frame_size)

# Creazione del video dai frame
for image_path in image_files:
    frame = cv2.imread(image_path)

    # Regola i livelli di luminosità e contrasto
    #alpha = 2
    #beta = 30
    #adjusted_frame = cv2.convertScaleAbs(frame, alpha=alpha, beta=beta)

    resized_frame = cv2.resize(frame, frame_size)
    out.write(resized_frame)

# Rilascio delle risorse
out.release()

