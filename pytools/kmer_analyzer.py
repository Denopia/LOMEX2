import sys
import math
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="darkgrid")

nuc_prio = {'C' : 0, 'A' : 1, 'T' : 2, 'G' : 3}
nuc_comp = {'C' : 'G', 'A' : 'T', 'T' : 'A', 'G' : 'C'}


def clean_line(line, counts=False):
	if line:
		if counts:
			return line.split()[0].strip(), 1
			return line.split()[0].strip(), int(math.ceil(float(line.split()[-1].strip())))
		else:
			return line.split()[0].strip()
	else:
		if counts:
			return "", 0
		else:
			return ""

def orient(kmer):
	if len(kmer) == 0:
		return ""
	complement = ""
	for nuc in kmer:
		if nuc not in nuc_prio:
			return ""
		complement = nuc_comp[nuc] + complement
	for i in range(len(kmer)):
		if nuc_prio[kmer[i]] == nuc_prio[complement[i]]:
			continue
		if nuc_prio[kmer[i]] > nuc_prio[complement[i]]:
			return kmer
		if nuc_prio[kmer[i]] < nuc_prio[complement[i]]:
			return complement
	return kmer


# THIS IS THE MOST RECENT AND EFFICIENT COMPARATOR. GIVE K-MERS IN ALPHABETICAL ORDER!
def compare(real_path, extracted_path, k):
	real = 0
	extracted = 0
	rex = 0
	lastreal = ""
	lastex = ""
	exdup = 0
	realdup = 0
	counter = 0
	extracted_count = 0
	count_correct = 100*[0]
	count_false = 100*[0]
	with open(real_path, 'r') as rfile, open(extracted_path, 'r') as xfile:
		curex, extracted_count = clean_line(xfile.readline(), True)
		curreal = clean_line(rfile.readline(), False)
		if len(curreal) == k:
			real += 1
		if len(curex) == k:
			extracted += 1
		while True:
			counter+=1
			if counter % 1000 == 0 and False:
				print("Real:",real, "Extracted:", extracted, "Shared:", rex)
				print(curreal)
				print(curex)
			if len(curreal) == 0 and len(curex) == 0:
				break
			# CASE 1
			if curreal == curex:
				rex += 1 # same k-mer
				lastreal = curreal # save last real
				lastex = curex # save last extracted
				
				curreal = clean_line(rfile.readline(), False) # read new real
				curex, extracted_count = clean_line(xfile.readline(), True) # read new extracted
				
				# examine real
				if  lastreal != curreal:
					if len(curreal) == k:
						real += 1
				else:
					realdup += 1
				# examine extracted
				if  lastex != curex:
					if len(curex) == k:
						extracted += 1
						if extracted_count < 100:
							count_correct[extracted_count] += 1							
				else:
					print(curex)
					exdup += 1
			# CASE 2
			elif (len(curreal) == k and curreal < curex) or len(curex) == 0:
				lastreal = curreal # save last real
				curreal = clean_line(rfile.readline(), False)
				# examine real
				if  lastreal != curreal:
					if len(curreal) == k:
						real += 1
				else:
					realdup += 1
			# CASE 3
			elif (len(curex) == k and curex < curreal) or len(curreal) == 0:
				lastex = curex # save last extracted
				curex, extracted_count = clean_line(xfile.readline(), True)
				if  lastex != curex:
					if len(curex) == k:
						extracted += 1
						if extracted_count < 100:
							count_false[extracted_count] += 1	
				else:
					print(curex)
					exdup += 1

	precision = rex/extracted
	recall = rex/real

	if precision+recall != 0:
		f1 = (2*precision*recall)/(precision+recall)
	else:
		f1 = 0

	print("#############################################################")
	print("File: ", extracted_path)
	print("Real k-mers:", real)
	print("Real k-mer duplicates:", realdup)
	print("Extracted k-mers:", extracted)
	print("Extracted k-mer duplicates:", exdup)
	print("Shared k-mers:", rex)
	print("Precision:", precision)
	print("Recall:",recall)
	print("F1 score:",f1)
	print("Correct counts:", count_correct)
	print("False counts:", count_false)
	print("#############################################################")

	#plt.title("LoMeX extracted k-mer stats")
	#plt.plot(count_correct[2:], linestyle="-", color="r", alpha=0.8, linewidth=2.0)
	#plt.plot(count_false[2:], linestyle="-", color="b", alpha=0.8, linewidth=2.0)
	#plt.legend(["Correct","False"])
	#plt.xlabel("Extracted k-mer count")
	#plt.ylabel("Extracted k-mer count count")
	#plt.show()


def read_next(rfile, xfile, read_real, read_extracted, real_kmer_count, extracted_kmer_count, extracted_duplicate_count, real_duplicate_count, curex_abundance, curreal_abundance, current_extracted, current_real):

	new_curex_abundance = curex_abundance 
	new_curreal_abundance = curreal_abundance
	
	new_extracted_duplicate_count = extracted_duplicate_count 
	new_real_duplicate_count = real_duplicate_count
		
	new_extracted_kmer_count = extracted_kmer_count
	new_real_kmer_count = real_kmer_count
	
	new_current_extracted = current_extracted 
	new_current_real = current_real
	
	if read_real:
		while True:
			next_real = clean_line(rfile.readline(), False)
			if next_real == "":
				new_current_real = ""
				break
			if next_real == current_real:
				new_real_duplicate_count += 1
			else:
				new_current_real = next_real
				new_real_kmer_count += 1
				break

	if read_extracted:
		while True:
			next_extracted, next_abundance = clean_line(xfile.readline(), True)
			if next_extracted == "":
				new_current_extracted = ""
				new_curex_abundance = 0
				break
			if next_extracted == current_extracted:
				new_extracted_duplicate_count += 1
			else:
				new_current_extracted = next_extracted
				new_curex_abundance = next_abundance
				new_extracted_kmer_count += 1
				break
	return new_current_extracted, new_current_real, new_curex_abundance, new_curreal_abundance, new_real_kmer_count, new_extracted_kmer_count, new_extracted_duplicate_count, new_real_duplicate_count


def compare_v2(real_path, extracted_path, k):
	extracted_kmer_count = 0
	real_kmer_count = 0
	real_extracted_kmer_count = 0

	extracted_duplicate_count = 0	
	real_duplicate_count = 0

	curex_abundance = 0
	curreal_abundance = 0

	current_extracted = "" 
	current_real = ""
	
	read_real = True
	read_extracted = True

	count_correct = 60*[0]
	count_false = 60*[0]

	counter = 0

	with open(real_path, 'r') as rfile, open(extracted_path, 'r') as xfile:
		while True:
			current_extracted, current_real, curex_abundance, curreal_abundance, real_kmer_count, extracted_kmer_count, extracted_duplicate_count, real_duplicate_count = read_next(rfile, xfile, read_real, read_extracted, real_kmer_count, extracted_kmer_count, extracted_duplicate_count, real_duplicate_count, curex_abundance, curreal_abundance, current_extracted, current_real)

			counter+=1
			#if (counter%100 == -1):
				#print("EXTRACTED:", current_extracted)
				#print("REAL:", current_real)
			# CASE 6
			if current_real == "" and current_extracted == "":
				break

			# CASE 1
			elif current_extracted == current_real:
				real_extracted_kmer_count += 1
				read_real = True
				read_extracted = True
				if curex_abundance < 60:
					count_correct[curex_abundance] += 1
								
			# CASE 2
			elif current_extracted < current_real and current_extracted != "":
				read_real = False
				read_extracted = True
				#print(current_extracted)
				#print(current_extracted, "EXTRACTED", curex_abundance)
				if curex_abundance < 60:
					count_false[curex_abundance] += 1
			
			# CASE 3
			elif current_extracted > current_real and current_real != "":
				#print(current_real, "REAL")				
				read_real = True
				read_extracted = False

			# CASE 4
			elif current_real == "" and current_extracted != "":
				read_real = False
				read_extracted = True
			
			# CASE 5
			elif current_extracted == "" and current_real != "":
				read_extracted = False
				read_real = True
			


	precision = real_extracted_kmer_count/extracted_kmer_count
	recall = real_extracted_kmer_count/real_kmer_count

	if precision+recall != 0:
		f1 = (2*precision*recall)/(precision+recall)
	else:
		f1 = 0

	print("#############################################################")
	print("First file (extracted k-mers): ", extracted_path)
	print("Second file (real k-mers): ", real_path)
	print("Real k-mers:", real_kmer_count)
	print("Real k-mer duplicates:", real_duplicate_count)
	print("Extracted k-mers:", extracted_kmer_count)
	print("Extracted k-mer duplicates:", extracted_duplicate_count)
	print("Correct extracted k-mers:", real_extracted_kmer_count)
	print("False extracted k-mers:", extracted_kmer_count - real_extracted_kmer_count)
	print("Missing real k-mers:", real_kmer_count - real_extracted_kmer_count)
	print("Precision:", precision)
	print("Recall:", recall)
	print("F1 score:", f1)
	print("Correct counts:", count_correct)
	print("False counts:", count_false)
	print("#############################################################")

	'''
	plt.title("DSK extracted k-mer stats, min abundance = 10")
	plt.plot(count_correct[2:], linestyle="-", color="r", alpha=0.8, linewidth=2.0)
	plt.plot(count_false[2:], linestyle="-", color="b", alpha=0.8, linewidth=2.0)
	plt.legend(["Correct","False"])
	plt.xlabel("Extracted k-mer count")
	plt.ylabel("Extracted k-mer count count")
	plt.show()
	'''


def compare_v3(real_path, extracted_path, k):
	extracted_kmer_count = 0
	real_kmer_count = 0
	real_extracted_kmer_count = 0

	extracted_duplicate_count = 0	
	real_duplicate_count = 0

	curex_abundance = 0
	curreal_abundance = 0

	current_extracted = "" 
	current_real = ""
	
	read_real = True
	read_extracted = True

	count_correct = 60*[0]
	count_false = 60*[0]

	counter = 0

	with open(real_path, 'r') as rfile, open(extracted_path, 'r') as xfile:
		while True:
			current_extracted, current_real, curex_abundance, curreal_abundance, real_kmer_count, extracted_kmer_count, extracted_duplicate_count, real_duplicate_count = read_next(rfile, xfile, read_real, read_extracted, real_kmer_count, extracted_kmer_count, extracted_duplicate_count, real_duplicate_count, curex_abundance, curreal_abundance, current_extracted, current_real)

			counter+=1
			if (counter%100 == -1):
				print("EXTRACTED:", current_extracted)
				print("REAL:", current_real)
			# CASE 6
			if current_real == "" and current_extracted == "":
				break

			# CASE 1
			elif current_extracted == current_real:
				real_extracted_kmer_count += 1
				read_real = True
				read_extracted = True
				if curex_abundance < 60:
					count_correct[curex_abundance] += 1
								
			# CASE 2
			elif current_extracted < current_real and current_extracted != "":
				read_real = False
				read_extracted = True
				if curex_abundance < 60:
					count_false[curex_abundance] += 1
			
			# CASE 3
			elif current_extracted > current_real and current_real != "":
				read_real = True
				read_extracted = False

			# CASE 4
			elif current_real == "" and current_extracted != "":
				read_real = False
				read_extracted = True
			
			# CASE 5
			elif current_extracted == "" and current_real != "":
				read_extracted = False
				read_real = True
			


	precision = real_extracted_kmer_count/extracted_kmer_count
	recall = real_extracted_kmer_count/real_kmer_count

	if precision+recall != 0:
		f1 = (2*precision*recall)/(precision+recall)
	else:
		f1 = 0

	print("#############################################################")
	print("First file (extracted k-mers): ", extracted_path)
	print("Second file (real k-mers): ", real_path)
	print("Real k-mers:", real_kmer_count)
	print("Real k-mer duplicates:", real_duplicate_count)
	print("Extracted k-mers:", extracted_kmer_count)
	print("Extracted k-mer duplicates:", extracted_duplicate_count)
	print("Correct extracted k-mers:", real_extracted_kmer_count)
	print("False extracted k-mers:", extracted_kmer_count - real_extracted_kmer_count)
	print("Missing real k-mers:", real_kmer_count - real_extracted_kmer_count)
	print("Precision:", precision)
	print("Recall:", recall)
	print("F1 score:", f1)
	print("Correct counts:", count_correct)
	print("False counts:", count_false)
	print("#############################################################")

	'''
	plt.title("DSK extracted k-mer stats, min abundance = 10")
	plt.plot(count_correct[2:], linestyle="-", color="r", alpha=0.8, linewidth=2.0)
	plt.plot(count_false[2:], linestyle="-", color="b", alpha=0.8, linewidth=2.0)
	plt.legend(["Correct","False"])
	plt.xlabel("Extracted k-mer count")
	plt.ylabel("Extracted k-mer count count")
	plt.show()
	'''

def main(real_path, extracted_path, k):
	print("Analyzing k-mers with")
	compare_v2(real_path, extracted_path, k)
	#compare(real_path, extracted_path, k)

if __name__ == "__main__":
   main(sys.argv[1], sys.argv[2], int(sys.argv[3]))

	
