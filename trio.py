
import vcf

def consistency(child_allele_count, mother_allele_count, father_allele_count):
	if child_allele_count == 0:
		return (mother_allele_count != 2 and father_allele_count != 2)
	elif child_allele_count == 1:
		return not (mother_allele_count == 2 and father_allele_count == 2) and not (mother_allele_count == 0 and father_allele_count == 0)
	else:
		return (mother_allele_count != 0) or (father_allele_count != 0)

def get_allele_count(gt, allele):
	return gt.count(str(allele))

def check_consistency (lina):
	consistent = []
	inconsistent = []
	child_gt = record.genotype('SAMP01')['GT'].split("/")
	mother_gt = record.genotype('SAMP02')['GT'].split("/")
	father_gt = record.genotype('SAMP03')['GT'].split("/")
	if len(record.ALT) > 1:
		for allele in xrange(len(record.ALT) + 1):
			child_allele_count = get_allele_count(child_gt, allele)
			mother_allele_count = get_allele_count(mother_gt, allele)
			father_allele_count = get_allele_count(father_gt, allele)
			if consistency(child_allele_count, mother_allele_count, father_allele_count):
				consistent.append(allele)
			else:
				inconsistent.append(allele)
	else:
		child_allele_count = get_allele_count(child_gt, 1)
		mother_allele_count = get_allele_count(mother_gt, 1)
		father_allele_count = get_allele_count(father_gt, 1)
		if consistency(child_allele_count, mother_allele_count, father_allele_count):
			consistent.append(1)
		else:
			inconsistent.append(1)
	return (consistent, inconsistent)


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

vcf_reader = vcf.Reader(open('test.vcf', 'r'))

for record in vcf_reader:

	for line in open('trio_skra.txt'):

		samples = line.rstrip().split("\t")
		child_gt = record.genotype(samples[0])['GT'].split("/")
		mother_gt = record.genotype(samples[1])['GT'].split("/")
		father_gt = record.genotype(samples[2])['GT'].split("/")

		print child_gt
		print mother_gt
		print father_gt

		#fallid check_consistency lista sem telur consistency. Thegar vid erum bara med einn reference tha er annadhvort 1 i consistent eda
		#1 i inconsistent. THhegar thad eru > 1 reference tha er thetta inconsistent ef thad er amk 1 inconsistent stak.
		(consistent, inconsistent) = check_consistency(record)

		if (father_gt[0] == father_gt[1]) and (mother_gt[0] == mother_gt[1]): #skoda hvort thetta se type1
			type_1 += 1
			#tjekkum hvort thad se consistent eda inconsistent
			if len(inconsistent) == 0:
				consistent_type1 += 1
			else:
				inconsistent_type1 += 1

		elif (father_gt[0] == father_gt[1]) or (mother_gt[0] == mother_gt[1]): #skoda hvort thetta se type2
			type_2 += 1
			#tjekkum hvort thad se consistent eda inconsistent
			if len(inconsistent) == 0:
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
			if len(inconsistent) == 0:
				consistent_type3 += 1
			else:
				inconsistent_type3 += 1
			if child_gt[0] == 0 and child_gt[1] == 0:
				o1 += 1
			if child_gt[0] == 1 and child_gt[1] == 1:
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

	print 'Variant info: {}'.format(record.ID)
	print 'type 1: {}'.format(type_1)
	print 'consistent type 1: {}'.format(consistent_type1)
	print 'inconsistent type 1: {}'.format(inconsistent_type1)
	print 'type 2: {}'.format(type_2)
	print 'consistent type 2: {}'.format(consistent_type2)
	print 'inconsistent type 2: {}'.format(inconsistent_type2)
	#tokum chi square test fyrir type 2 ef thad eru einhverjir type 3 trio
	if type_2 > 0:
		chi_statistic_type2 = (heterozygous_type2-type_2*0.5)**2/(type_2*0.5) + (homozygous_type2-type_2*0.5)**2/(type_2*0.5)
		print 'chi square statistic fyrir type 2: {}'.format(chi_statistic_type2)
	print 'type 3: {}'.format(type_3)
	print 'consistent type 3: {}'.format(consistent_type3)
	print 'inconsistent type 3: {}'.format(inconsistent_type3)
	#tokum chi square test fyrir type 3 ef thad eru einhverjir type 3 trios
	if type_3 > 0:
		chi_statistic_type3 = (o1-type_3*0.25)**2/(type_3*0.25) + (o2-type_3*0.5)**2/(type_3*0.5) + (o3-type_3*0.25)**2/(type_3*0.25)
		print 'chi square statistic fyrir type 3: {}'.format(chi_statistic_type3)
	print '\n'
	print '\n'

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
	




