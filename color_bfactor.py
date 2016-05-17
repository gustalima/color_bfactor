from pymol import cmd
import os

def colorbfactor(*args):
	"""
	DESCRIPTION

	    Color objects based on their sequence alignment. User must install clustalo (from repository,
	    if available) and the zero_residues plugin http://www.pymolwiki.org/index.php/Zero_residues

	USAGE
	    Select the objects you want to work in PyMOL's object list and then:

	    colorbeta
        """


	#Start with base color and get selected objects in PyMOL
	cmd.do('color white')
	enabled = cmd.get_names(enabled_only=1)

	#Renumber all residues using zero_residues module for pymol http://www.pymolwiki.org/index.php/Zero_residues
	for i in enabled:
		cmd.do('zero_residues %s, 1'%i)

	#Save fasta sequences for alignment
	cmd.do('save sequences.fasta, %s' %' '.join(enabled))

	#Create alignment using Clustal O
	os.system('clustalo -i sequences.fasta -outfmt=selex --wrap=10000')
	entry = open('clustal.aln','r').readlines()[3:]

	os.system('clustalo -i sequences.fasta -outfmt=selex --wrap=70')
	fixed_entry = open('clustal.aln','r').readlines()[3:]

	#Get the first entry of alignment itself excluding sequence identifier
	start = entry[0].index(entry[0].split()[1])

	#Create sequences and id lists removing empty spaces
	tags = []
	seqs = []
	for i in entry:
		tags.append(i[:start].replace(' ',''))
		seqs.append(i[start:])

	#Align objects in PyMOL
	cmd.do('alignto %s; center' %tags[0])

	#Remove empty entries in list and remove \n from end of lines
	tags = [x for x in tags if len(x)!=0]
	cons = seqs[-1].replace('\n','')

	seqs = [x for x in seqs if len(x)!=0]
	seqs = [x.replace('\n','') for x in seqs[:-1]]


	#Print a very basic checkup for problems with sequences and ID. Makes no checks on sequence content
	if len(tags) == len(seqs):
		print 'Pre-processing OK. Now editing PyMOL objects\n'
	else:
		print 'Ops.. Error with your sequences. Try to change sequences names\n'

	#Print the color and symbol pattern for output.
	print '( ) for non-conservative AA  -> Red'
	print '(.) for semi-conservative AA -> Orange'
	print '(:) for conservative AA      -> Yellow'
	print '(*) for identical AA         -> Light grey'
	print '\n'

	for i in fixed_entry:
		print i.replace('\n','')

	#Create sequence lists starting from the first AA instead of insertion, if applies
	seqs_blunt = []
	for i in seqs:
		helper = i
		while helper[0]=='-':
			helper = helper[1:]
		while helper[-1]=='-':
			helper = helper[:-1]
		seqs_blunt.append(helper)

	#Matches the conservation to sequences scheme
	cons_blunt = []
	for i,j in zip(seqs_blunt,seqs):
		cons_blunt.append( cons[j.index(i):])

	#Solve the problem of insertion in middle of sequences
	conservation = []
	for i,j in zip(seqs_blunt,cons_blunt):
		counter = 0
		helper = ''
		for k,l in zip(i,j):

			if k=='-':
				helper +=''
			else:
				helper += l
			counter +=1
		conservation.append(helper)

	#Recolor selected objects in PyMOL
	for i,j in zip(tags, conservation):
		for k in enumerate(j):
			if k[1] == '*':
				cmd.color('grey80', "resi %s and %s" %(str(k[0]+1),i))
			if k[1] == ':':
				cmd.color('olive' , "resi %s and %s" %(str(k[0]+1),i))
			if k[1] == '.':
				cmd.color('orange', "resi %s and %s" %(str(k[0]+1),i))
			if k[1] == ' ':
				cmd.color('red'   , "resi %s and %s" %(str(k[0]+1),i))

	####### RMSD calculation #########


	if len(enabled)==2:

		cmd.show_as("cartoon")
		cmd.cartoon("putty")
		print "Calculating RMSD for sidechains. It ONLY works for two objects."
		print "Make sure you don't have more than two selected"
		import numpy
		def dist(x,y):
			#rmsd of two residues
			x = numpy.array(x)
			y = numpy.array(y)
			return numpy.sqrt(numpy.sum((x-y)**2))

		def centroid(lista):
			#return the centroid of the sidechain for a given xyz list
			alista = numpy.array(lista)
			length = alista.shape[0]
			sum_x = numpy.sum(alista[:, 0])
			sum_y = numpy.sum(alista[:, 1])
			sum_z = numpy.sum(alista[:, 1])
			return sum_x/length, sum_y/length, sum_z/length


		xyz_table = []
		for i,j in zip(tags,seqs):
			helper = []
			for k in enumerate(j.replace('-','')):

				helper.append(cmd.get_model('resi %s and %s'%(str(k[0]+1),i),1).get_coord_list())

			xyz_table.append(helper)

		centroid_table = []
		print 'Calculating centeroids'
		for i in xyz_table:
			helper = []
			for j in i:
				helper.append(list(centroid(j)))
			centroid_table.append(helper)
		dist_table = []

		print 'Done!'

		first 		 = ''
		second		 = ''
		centr_first  = ''
		centr_second = ''
		for i,j in zip(range(len(seqs[0])),range(len(centroid_table[0]))):
			if seqs[0][i]=='-' or seqs[1][i]=='-':
				pass
			else:
				first   += seqs[0][i]
				second  += seqs[1][i]
		equal_seqs 		= [first,second]
		equal_centroids = [centr_first,centr_second]
		ini	= seqs[1].index(equal_seqs[1][:4])
		fin	= seqs[1][::-1].index(equal_seqs[1][::-1][:4])
		if len(centroid_table[1])==len(centroid_table[0]):
			pass
		elif len(centroid_table[1])>len(centroid_table[0]):
			centroid_table[1] = centroid_table[1][ini:-fin]
		else:
			centroid_table[0] = centroid_table[0][ini:-fin]
		rmsd = []
		print 'Calculating RMSDs'
		for i,j in zip(centroid_table[0],centroid_table[1]):
			rmsd.append(dist(i,j))
		print 'Done!'
		rmsd = [10*x/max(rmsd) for x in rmsd]
		diff =  len(seqs[0].replace('-','')) - len(seqs[1].replace('-',''))
		addone = 0
		addtwo = 0

		if diff > 0:
			addtwo = diff
		elif diff < 0:
			addone = -diff
		else:
			pass
		for j in enumerate(rmsd):
			cmd.alter('%s and resi %s and n. CA'%(tags[1],str(j[0]+1+addone)), 'b=%s'%j[1])
			cmd.alter('%s and resi %s and n. CA'%(tags[0],str(j[0]+1+addtwo)), 'b=%s'%j[1])









cmd.extend("color_putty", colorbfactor);
colorbfactor()
