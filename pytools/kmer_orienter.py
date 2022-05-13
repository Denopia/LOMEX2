import sys

nuc_prio = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3}
nuc_comp = {'C' : 'G', 'A' : 'T', 'T' : 'A', 'G' : 'C'}


def orient_file(old_path, new_path, k, a):
	#l = 0
	with open(old_path, 'r') as rfile, open(new_path, 'w') as wfile:
		for line in rfile:
			if int(line.split()[-1]) < a:
				continue
			oriented_line = orient(clean_line(line))
			if len(oriented_line) != k:
				print("Mystery length error\n")
				continue
			wfile.write(oriented_line+"\n")
	
def clean_line(line):
	if line:
		return line.split()[0].strip()
	else:
		return ""

def orient(kmer):
	complement = ""
	for nuc in kmer:
		if nuc not in nuc_prio:
			return "Z"
		complement = nuc_comp[nuc] + complement
	for i in range(len(kmer)):
		if nuc_prio[kmer[i]] == nuc_prio[complement[i]]:
			continue
		if nuc_prio[kmer[i]] < nuc_prio[complement[i]]:
			return kmer
		if nuc_prio[kmer[i]] > nuc_prio[complement[i]]:
			return complement
	return kmer


def main(old_path, new_path, k, a):
	print("Orienting k-mer file\n")
	orient_file(old_path, new_path, k, a)

if __name__ == "__main__":
   main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))


