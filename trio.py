from __future__ import division
from __future__ import print_function
import vcf
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("vcf", type = str)
parser.add_argument("--trio", type = str)
parser.add_argument("--twins", type = str)
parser.add_argument("--info_list", type = str, nargs = '*')

args = parser.parse_args()

def consistency(child_allele_count, mother_allele_count, father_allele_count):
	if child_allele_count == 0:
		return (mother_allele_count != 2 and father_allele_count != 2)
	elif child_allele_count == 1:
		return not (mother_allele_count == 2 and father_allele_count == 2) and not (mother_allele_count == 0 and father_allele_count == 0)
	else:
		return (mother_allele_count != 0) or (father_allele_count != 0)

def get_allele_count(gt, allele):
	return gt.count(str(allele))

def check_consistency (record, samples, target_allele):
	consistent = []
	inconsistent = []
	child_gt = record.genotype(samples[0])['GT'].split("/")
	mother_gt = record.genotype(samples[1])['GT'].split("/")
	father_gt = record.genotype(samples[2])['GT'].split("/")

	child_allele_count = get_allele_count(child_gt, target_allele)
	mother_allele_count = get_allele_count(mother_gt, target_allele)
	father_allele_count = get_allele_count(father_gt, target_allele)
	if consistency(child_allele_count, mother_allele_count, father_allele_count):
		return True
	else:
		return False

def statistic_tests(record, fjoldi_samples, target_allele, items):

	type_1 = 0
	consistent_type1 = 0
	inconsistent_type1 = 0
	heterozygous_type1 = 0
	homozygous_type1 = 0

	type_2 = 0
	consistent_type2 = 0
	inconsistent_type2 = 0
	heterozygous_type2 = 0
	homozygous_type2 = 0

	type_3 = 0
	consistent_type3 = 0
	inconsistent_type3 = 0
	o1 = 0
	o2 = 0
	o3 = 0

	P = 0
	H = 0
	Q = 0
	harwey_chi = 0


	for sample in vcf_reader.samples:
		allele =  record.genotype(sample)['GT'].split("/")

		if "." in allele:
			continue

		count = allele.count(str(target_allele))
		if count == 0:
			P += 1
		elif count == 1:
			H += 1
		else:
			Q += 1

	fjoldi_samples = P + H + Q
	p = (P + H/2)/fjoldi_samples
	q = (Q + H/2)/fjoldi_samples

	expected_1 = (p * p)*fjoldi_samples #expected fjoldi af 0/0
	expected_2 = (2*q*p)*fjoldi_samples #expected fjoldi af 1/0
	expected_3 = (q*q)*fjoldi_samples   #expected fjoldi af 1/1

	if expected_1 > 0 and expected_2 > 0 and expected_3 > 0:
		harwey_chi = (P-expected_1)**2/expected_1 + (H-expected_2)**2/expected_2 + (Q-expected_3)**2/expected_3

	o1 = 0
	o2 = 0
	o3 = 0

	if args.trio is not None:
		for line in open(args.trio):

			samples = line.rstrip().split("\t")
			child_gt = record.genotype(samples[0])['GT'].split("/")
			mother_gt = record.genotype(samples[1])['GT'].split("/")
			father_gt = record.genotype(samples[2])['GT'].split("/")

			if "." in child_gt or "." in mother_gt or "." in father_gt:
				continue

			#fallid check_consistency lista sem telur consistency. Thegar vid erum bara med einn reference tha er annadhvort 1 i consistent eda
			#1 i inconsistent. THhegar thad eru > 1 reference tha er thetta inconsistent ef thad er amk 1 inconsistent stak.
			consistency = check_consistency(record, samples, target_allele)

			f1 = father_gt[0] == str(target_allele)
			f2 = father_gt[1] == str(target_allele)
			m1 = mother_gt[0] == str(target_allele)
			m2 = mother_gt[1] == str(target_allele)

			if (f1 == f2) and (m1 == m2): #skoda hvort thetta se type1
				type_1 += 1
				#tjekkum hvort thad se consistent eda inconsistent
				if consistency:
					consistent_type1 += 1
				else:
					inconsistent_type1 += 1

			elif (f1 == f2) or (m1 == m2): #skoda hvort thetta se type2
				type_2 += 1
				#tjekkum hvort thad se consistent eda inconsistent
				if consistency:
					consistent_type2 += 1
					#skodum hversu morg type2 eru heterozygous og hversu morg eru homozygous. ATH inconsistent type 2 flokkad sem heterozygous.
					if child_gt[0] == child_gt[1]:
						homozygous_type2 += 1
					else:
						heterozygous_type2 +=1
				else:
					inconsistent_type2 += 1
					heterozygous_type2 +=1

			else: #annars er thetta type 3
				type_3 += 1
				if consistency:
					consistent_type3 += 1
				else:
					inconsistent_type3 += 1
				if child_gt[0] != str(target_allele) and child_gt[1] != str(target_allele):
					o1 += 1
				elif child_gt[0] == str(target_allele) and child_gt[1] == str(target_allele):
					o3 += 1
				else:
					o2 +=1

	#info um chi square tilgatuprofid:
	#i badum tilvikum er nulltilgatan: No significant difference between observed and expected og gagntilgatan er
	# there is a significant difference between observed and expected svo ef vid neitum nulltilgatunni tha getum vid sagt
	# med X mikilli vissu (eftir thvi hvad oryggisbilid er) ad gognin seu ekki rett dreyfd.
	#fyrir type 2 tha skoda eg hversu morg born eru heterozygous og hversu morg born eru homozygous. Inconsistent
	#flokkast sem heterozygous. Expected er n*0,5 hja badum og svo tel eg hversu margir eru homo og hversu margir
	#eru heter og framkvaemi svo nulltilgatuprofid.
	#fyrir type 3 tha set eg o1 = 0/0 og expected n*0.25, o2 = 1/0 og expected n*0,5 og o3 = 1/1 og expected n*0.25
	#framkvaemi svo tilgatuprofid.

	if type_2 > 0:
		chi_statistic_type2 = (heterozygous_type2-type_2*0.5)**2/(type_2*0.5) + (homozygous_type2-type_2*0.5)**2/(type_2*0.5)
	else:
		chi_statistic_type2 = 'n/a' #thvi vid erum tha ad deila med nulli, spurning um tulkun a thvi 

	if type_3 > 0:
		chi_statistic_type3 = (o1-type_3*0.25)**2/(type_3*0.25) + (o2-type_3*0.5)**2/(type_3*0.5) + (o3-type_3*0.25)**2/(type_3*0.25)
		
	else:
		chi_statistic_type3 = 'n/a' #thvi vid erum tha ad deila med nulli, spurning um tulkun a thvi 


	if args.twins is not None:
		twins_passar = 0
		twins_passar_ekki = 0
		twins_prosenta = 0
		for lina in open(args.twins):
			samples_twins = lina.rstrip().split("\t")
			twin1 = record.genotype(samples_twins[0])['GT'].split("/")
			twin2 = record.genotype(samples_twins[1])['GT'].split("/")
			key = samples_twins[0] + samples_twins[1]

			if key not in mz_twins:
				mz_twins[key] = [0, 0]

			if twin1.count('1') == twin2.count('1'):
				twins_passar += 1
				mz_twins[key][0] += 1
			else:
				twins_passar_ekki += 1
				mz_twins[key][1] += 1

		
		twins_prosenta = twins_passar/(twins_passar + twins_passar_ekki)

		file.write(record.ID + ':' + str(target_allele) + '\t' + str(twins_prosenta))
		file.write('\n')

	gildi = []
	if args.info_list is not None:
		for bla in items:
			if bla in record.INFO:
				if type(record.INFO[bla]) is list:
					if len(record.INFO[bla]) == len(record.ALT):
						gildi.append(record.INFO[bla][target_allele-1])
					elif len(record.INFO[bla]) == len(record.ALT) + 1:
						gildi.append(record.INFO[bla][target_allele])
					else:
						gildi.append(record.INFO[bla])
				else: 
					gildi.append(record.INFO[bla])
			else:
				gildi.append('n/a')


	#hendum ollum upplysingunum inn i global dictionary sem geymir thaer
	trio_listi = []
	info_listi = []
	if args.trio is not None:
		trio_listi = [type_1, consistent_type1, inconsistent_type1, type_2, consistent_type2, inconsistent_type2, chi_statistic_type2, type_3, consistent_type3, inconsistent_type3, chi_statistic_type3]

	prent_var = trio_listi + gildi

	print (record.ID + ':' + str(target_allele), harwey_chi, *prent_var, sep = "\t" )
#-----------------------------------------------------------------------------------------------------------------


vcf_reader = vcf.Reader(open(args.vcf, 'r'))

mz_twins = {} # Map with first monozygote twin as key, and value is a list where the first value are number of constistencies
 # and the second is the number of inconsistensies

variants = {}

items = []
prent_trio = []

if args.info_list is not None:
	items = args.info_list

if args.trio is not None:
	prent_trio = ['#Type 1', "# Consistent Type1", "Inconsistent Type1", "# Type 2", "# Consistent Type2", "Inconsistent Type2", "Type 2 Chi Square Statistic" "# Type 3", "# Consistent Type3", "Inconsistent Type3" "Type 3 Chi Square Statistic"]

prenti_listi = prent_trio + items

print ('ID', "Variant H.W C.S.S", *prenti_listi, sep = "\t" )

fjoldi_samples = len(vcf_reader.samples)

snippur = 0
deletion = 0
insertion = 0
complx = 0

if args.twins is not None:
	file = open ("twin_output.txt", "w")
	file.write('ID' + '\t' + 'Hlutfall')
	file.write ('\n')


for record in vcf_reader:

	if len(record.ALT) == 1:
		statistic_tests(record, fjoldi_samples, 1, items)
	else:
		for i in range(len(record.ALT)):
			statistic_tests(record, fjoldi_samples, i+1, items)

	allele_len_1 = sum(len(i) == 1 for i in record.ALT) + int(len(record.REF) == 1)

	if allele_len_1 == (len(record.ALT) + 1):
		snippur += 1
	elif allele_len_1 != 1:
		complx += 1
	elif len(record.REF) == 1:
		insertion += 1
	else:
		deletion +=1

if args.twins is not None:
	file.close()

print ('Snippur', 'Insertion', 'Deletion', 'Complex', sep = '\t')
print (snippur, insertion, deletion, complx, sep = '\t')




	




