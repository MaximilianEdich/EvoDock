# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 14:51:42 2017

Give it a textfile which contains commands for the command line and it
will run trough each line and will do the command automatically.

@author: Maximilian Edich
"""

import sys
import os
import time

"""
Takes a time in seconds. Prints it depending on the size. Times less than
5sec are displayed in millis. Times up to a minute are displayed just in
seconds. All bigger times are displayed in minutes and seconds.
"""


def printTime(timeSec):
    # if the time span is big, display it in seconds and minutes, otherwise in millis
    if (timeSec >= 5):
        if (timeSec >= 60):
            m, s = divmod(timeSec, 60)
            print("Runtime: " + str(int(m)) + " min " + str(round(s)) + " sec")
        else:
            print("Runtime in seconds: " + str(round(timeSec)))
    else:
        print("Runtime in millis: " + str(round(timeSec * 1000)))


if len(sys.argv) >= 2:
    if os.path.isfile(sys.argv[1]):
        start = time.time()

        file = open(sys.argv[1], 'r')
        lines = file.readlines()
        file.close()

        # calculate number of real commands and ignore whitespaces
        numOfCommands = 0
        for line in lines:
            if line != '\n':
                numOfCommands += 1

        # run commands
        command = 0
        for line in lines:
            command += 1
            print("Command " + str(command) + "/" + str(numOfCommands), end='\r')
            print("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - \n")
            print(line)
            try:
                os.system(line)
            except KeyboardInterrupt:
                exit("KEYBOARD INTERRUPT")
            except Exception:
                exit("ERROR")
        print("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        end = time.time()
        printTime(end - start)
        print("commands done")

    else:
        sys.exit("ERROR: Input has to be folder!")
else:
    sys.exit("ERROR: No Input!")
