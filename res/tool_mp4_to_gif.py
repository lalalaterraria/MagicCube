from moviepy.editor import *

v = VideoFileClip("res/demo.mp4")
v = v.resize((v.size[0]//2,v.size[1]//2))
v.subclip(0, 10).write_gif("res/1.gif")
v.subclip(10, 20).write_gif("res/2.gif")
v.subclip(20, 30).write_gif("res/3.gif")
v.subclip(147, 157).write_gif("res/4.gif")