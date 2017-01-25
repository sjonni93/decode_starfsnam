

import vcf

snippur = 0
deletion = 0
insertion = 0
complx = 0
villa = 0

vcf_reader = vcf.Reader(open('test.vcf', 'r'))
for record in vcf_reader:
	if len(record.ALT) > 1:
		continue
	elif len(record.REF) == 1 and len(record.ALT[0]) == 1:
		snippur += 1
	elif len(record.REF) == 1 and len(record.ALT[0]) > 1:
		insertion += 1
	elif len(record.REF) > 1 and len(record.ALT[0]) == 1:
		deletion += 1
	elif len(record.REF) > 1 and LEN(RECORD.alt[0]) == 1:
		complex += 1
	else:
		villa += 1

print snippur
print insertion
print deletion
print complx
print villa





