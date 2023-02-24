import random
import sys

n=int(sys.argv[1])

with open(sys.argv[2], "w") as f:
    for i in range(n):
        f.write(f"{random.random()*10000},{random.random()*10000},{random.random()*10000},{random.random()*100},{random.random()*100},{random.random()*100},{random.random()*1000}\n")

