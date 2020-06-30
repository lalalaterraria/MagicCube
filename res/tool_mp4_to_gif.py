from moviepy.editor import *

v = VideoFileClip("res/demo.mp4")
VideoFileClip("res/demo.mp4").subclip(0, 10).write_gif("res/1.gif")
VideoFileClip("res/demo.mp4").subclip(10, 20).write_gif("res/2.gif")
VideoFileClip("res/demo.mp4").subclip(20, 30).write_gif("res/3.gif")
VideoFileClip("res/demo.mp4").subclip(147, 157).write_gif("res/4.gif")