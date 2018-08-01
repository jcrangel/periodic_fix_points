#This programa helps to covert equations from mathematica (as plain text)
#to a maneageable form to copy into Matlab or other language.
#How to: Copy the mathematica cell as plaintext and paste it in file, use 
#this file as input of this script. 
import sys, getopt, re


# Get the parameters
def main(argv):
   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print('test.py -i <inputfile> -o <outputfile>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('test.py -i <inputfile> -o <outputfile>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
###############

   #regex = r"([a-zA-Z]+) (\d+)"
   fileIn=open(inputfile,"r")
   stringg=fileIn.read()
   regext="(\[t\])"
   #Weirly, for match one \ its needed four \\\\
   regexBrackets="\\\\\[(.*?)\]" #quita los \[Algo]
   regexSpaceMult="([a-z|A-Z])([\s]+)([a-z|A-Z])" #pone * donde hay multiplicaciones

   reg=re.sub(regext, r" ", stringg)
   reg=re.sub(regexBrackets, r"\1", reg)
   reg=re.sub(regexSpaceMult, r"\1*\3", reg)
   reg=re.sub(regexSpaceMult, r"\1*\3", reg)
   fileOut=open(outputfile,"w")
   fileOut.write(reg)


if __name__ == "__main__":
   main(sys.argv[1:])