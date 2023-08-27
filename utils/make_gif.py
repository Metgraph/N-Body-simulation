from PIL import Image

image_paths = [f'images/2022-{i}.png' for i in range(1000)]

# Apri le immagini
images = [Image.open(path) for path in image_paths]
output_path = "output.gif"

# Salva la GIF
print(f"Writing the gif to: {output_path}")

images[0].save(output_path, save_all=True, append_images=images[1:], duration=33, loop=0)

print("Completed!")

