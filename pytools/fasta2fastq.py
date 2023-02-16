import sys

def main(path):

        with open(path, "r") as rfile, open(path+".fastq", 'w') as wfile:
                fresh = True
                last_title = "nessu"
                readlen = 0
                for line in rfile:
                        line = line.strip()
                        if len(line) == 0:
                                continue
                        if line[0] == ">":
                                if not fresh:
                                        wfile.write("\n")
                                        wfile.write("+"+last_title[1:]+"\n")
                                        wfile.write(readlen*":"+"\n")
                                wfile.write("@"+line[1:]+"\n")
                                last_title = line
                                readlen = 0
                                fresh = False
                        else:
                             	wfile.write(line)
                                readlen += len(line)
                wfile.write("\n")
                wfile.write("+"+last_title[1:]+"\n")
                wfile.write(readlen*":"+"\n")

if __name__ == "__main__":
        main(sys.argv[1])








