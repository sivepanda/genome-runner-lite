import os

feature = "a"
genome = "b"

print(os.path.join(os.getcwd(),"track", feature + "_" + genome + ".bed")) 
