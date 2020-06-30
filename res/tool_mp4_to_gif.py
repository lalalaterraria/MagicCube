from moviepy.editor import *
begin, end = 20, 30
output = "res/3.gif"
VideoFileClip("res/demo.mp4").subclip(begin, end).write_gif(output)