import sys

def clean_seqline(seqline):
	seqline_parts = seqline.split()
	upper_seq = seqline_parts[0]
	no_deletions = upper_seq.replace('-', '')
	return no_deletions



def split(seq_file, new_file, k):
	with open(seq_file, 'r') as rfile, open(new_file, 'w') as wfile:
		for seqline in rfile.readlines():
			cleaned_seq = clean_seqline(seqline)
			if len(cleaned_seq) < k:
				continue
			i = 0
			while i + k <= len(cleaned_seq):
				wfile.write(cleaned_seq[i:i+k].upper()+" 999\n")				
				i+=1


def main(seq_file, new_file, k):
	split(seq_file, new_file, int(k))


if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2], sys.argv[3])
			
