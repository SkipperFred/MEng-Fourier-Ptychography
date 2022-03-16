from neopixel import Neopixel
from utime import sleep
import sys




numpx = 64
pixels = Neopixel(numpx,0,28,'GRB')

pixels.brightness(255)

pixels.set_pixel(0,(0,255,0))
pixels.show()
sleep(1)
pixels.set_pixel(0,(0,0,0))
pixels.show()

def LEDOn(LEDnum):
    pixels.set_pixel(LEDnum,(0,255,0))
    pixels.show()
    
def LEDOff(LEDnum):
    pixels.set_pixel(LEDnum,(0,0,0))
    pixels.show()
    
def LEDBright(brightnum):
    pixels.brightness(brightNum)


"""while True:
    
    data = sys.stdin.readline()
    
    #pippo = data.decode("utf-8")

    if data.len() > 0:
        pixels.set_pixel(32,(0,255,0))
        pixels.show()

    else:
        pixels.set_pixel(32,(0,255,0))
        pixels.show()"""