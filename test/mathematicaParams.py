#This programa helps to covert the asociations of paremeters
#from mathematica to matlab 

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
   regext="(->)"
   regexBrComments="\(\*[\s\S]*?\*\)"
   regexBrackets="\\\\\[(.*?)\]" #quita los \[Algo]
   regexComa=","
   regexvars="([a-zA-Z]+?)[\s]*\=" #save the variables for printing

   reg=re.sub(regext, r"=", stringg)
   reg=re.sub(regexBrComments, r"", reg)
   reg=re.sub(regexBrackets, r"\1", reg)
   reg=re.sub(regexComa, r";", reg)

   variables=re.findall(regexvars,reg)
   vars=""
   for variable in variables:
      vars=vars+" "+variable

   print(vars)   
   fileOut=open(outputfile,"w")
   fileOut.write(reg)


if __name__ == "__main__":
   main(sys.argv[1:])