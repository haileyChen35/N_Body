# pip install pdf2image pillow
# brew install poppler


from pdf2image import convert_from_path
from PIL import Image

# Convert PDF pages to images
pdf_path = "./solar.pdf"  # Change this to the correct file path
images = convert_from_path(pdf_path)

# Save as an animated GIF
gif_path = "solar_animation.gif"
images[0].save(gif_path, save_all=True, append_images=images[1:], duration=150, loop=0)

print(f"GIF saved as {gif_path}")
