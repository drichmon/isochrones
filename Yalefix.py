import os
import sys
import fileinput


'''the lines in yale files that contain the ages cause problems with numpy in Yalepreparer.py, because they are not commented like they are in other models. This script will comment those lines out.'''




path = '/home/richmond/V2/interpolated/' #location of your yale files

comment = "#"

term = "age(Gyr)="



'''
#comment age lines
for file in os.listdir(path): #puts a '#' in front of the ages
    file1 = path + '%s' % (file)
    for line in fileinput.input(file1, inplace = True):
        if term in line:
            print comment + "%s" % (line),
        else:
            print line,
         
'''



#reverse 
'''
for file in os.listdir(path): #takes # out of lines
    file1 = path + '%s' % (file)
    for line in fileinput.input(file1, inplace = True):
        if term in line:
            x = line.replace(comment, '')
            print x, 
        else:
            print line,

'''
#puts a '#' in front of the first line
'''
for file in os.listdir(path): 
    file1 = path + '%s' % (file)
    i = 1
    for line in fileinput.input(file1, inplace = True):
        if i < 2:
            if "#" in line:
                print line,
                i += 1 
            elif '#' not in line:
                print comment + "%s" % (line),
                i += 1
        else:
            print line, 
'''
'''
#places a '#' on the second line of the file
for file in os.listdir(path):
    file1 = path + '%s' % file
    i = 1
    for line in fileinput.input(file1, inplace = True):
        if i ==2:
            print comment + '%s' % line,
            i += 1
        else:
            print line,
            i += 1
'''


#comment one of the for loops if you need to only perform one action
