
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
type_2 = 0
consistent_type2 = 0
inconsistent_type2 = 0
type_3 = 0
consistent_type3 = 0
inconsistent_type3 = 0

vcf_reader = vcf.Reader(open('test.vcf', 'r'))

for record in vcf_reader:

	child_gt = record.genotype('SAMP01')['GT'].split("/")
	mother_gt = record.genotype('SAMP02')['GT'].split("/")
	father_gt = record.genotype('SAMP03')['GT'].split("/")

	print child_gt
	print mother_gt
	print father_gt

	(consistent, inconsistent) = check_consistency(record)

	if (father_gt[0] == father_gt[1]) and (mother_gt[0] == mother_gt[1]):
		type_1 += 1
		if len(inconsistent) == 0:
			consistent_type1 += 1
			print 'consistent_type1'
		else:
			inconsistent_type1 += 1
			print 'inconsistent_type1'
	elif (father_gt[0] == father_gt[1]) or (mother_gt[0] == mother_gt[1]):
		type_2 += 1
		if len(inconsistent) == 0:
			consistent_type2 += 1
			print 'consistent_type2'
		else:
			inconsistent_type2 += 1
			print 'inconsistent_type2'
	else:
		type_3 += 1
		if len(inconsistent) == 0:
			consistent_type3 += 1
			print 'consistent_type2'
		else:
			inconsistent_type3 += 1
			print 'inconsistent_type2'


print 'type 1: {}'.format(type_1)
print 'consistent type 1: {}'.format(consistent_type1)
print 'inconsistent type 1: {}'.format(inconsistent_type1)
print 'type 2: {}'.format(type_2)
print 'consistent type 2: {}'.format(consistent_type2)
print 'inconsistent type 2: {}'.format(inconsistent_type2)
print 'type 3: {}'.format(type_3)
print 'consistent type 3: {}'.format(consistent_type3)
print 'inconsistent type 3: {}'.format(inconsistent_type3)
	




