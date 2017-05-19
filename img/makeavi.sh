ffmpeg  -f image2 -framerate 24 -pattern_type sequence -start_number 0  -r 3 -i image%d.ppm -s 720x480 test.avi
