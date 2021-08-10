import sys

def clean_file(old_path, new_path, k, a):
	with open(old_path, 'r') as rfile, open(new_path, 'w') as wfile:
		for line in rfile:
			if int(line.split()[-1]) < a:
				continue
			
			if len(line.split()[0]) != k:
				continue
			wfile.write(line.split()[0]+"\n")
	
def clean_line(line):
	if line:
		return line.split()[0].strip()
	else:
		return ""

def main(old_path, new_path, k, a):
	print("Cleaning k-mer file\n")
	clean_file(old_path, new_path, k, a)
	#orient_file2(old_path, new_path, k)

if __name__ == "__main__":
   main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))


