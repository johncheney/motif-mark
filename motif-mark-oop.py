#! /usr/bin/env python 

import argparse
import re
import cairo

def get_args():             
    """This function passes arguments to a python script to identify binding motifs visually"""
    parser = argparse.ArgumentParser(description="A program to identify binding motifs visually")
    parser.add_argument("-f", "--filename", help="Your filename", required=True)               
    parser.add_argument("-m", "--motif", help="Your motif filename", required=True)               

    return parser.parse_args()

args = get_args()    
filename = (args.filename)
motif = (args.motif)


color_pal_list = [(1.0000, 0.0000, 0.0000),
                  (1.0000, 0.5478, 0.0000),
                  (0.4173, 1.0000, 0.0000),
                  (0.0000, 0.0174, 1.0000),
                  (0.3034, 0.0000, 1.0000),
                  (1.0000, 1.0000, 0.0000)]

fasta_dict={}
fasta_list=[]

pattern = "[a-z]+"
patternUPPER = "[A-Z]+"

gene_dict={}

gene_objects = {}
class Gene:
    def __init__(self, header, preintron, exon, postintron, total_length, full_gene_list):
            '''The Gene class'''

            self.header = header
            self.preintron = preintron
            self.exon = exon
            self.postintron = postintron
            self.total_length = total_length
            self.full_gene_list = full_gene_list
            self.motif_position_dict = {}

temp_file = "./temp.fasta"

with open(temp_file, "w") as temp:
    with open(filename, "r") as fh:
        firstline = True			            #opening the fastq file to read it 
        for line in fh:                         
            if line[0] == '>':                  
                header = line.strip() #removing end of line whitespace characters
                if firstline == False:
                    temp.write('\n')
                temp.write(header)
                temp.write('\n')
                firstline = False
            else:
                seq = line.strip()
                temp.write(seq)
        temp.write('\n')
        temp.close() 
    

def genie(): 
    with open(temp_file , "r") as gfh:
        for item in gfh:
            Gene.header = item.strip()
            item = gfh.readline()
            findall = re.findall(pattern,item) # List of Match objects.
            findall_upper= re.findall(patternUPPER,item) # List of Match objects.
            introns_list = findall
            exon_list = findall_upper
            Gene.full_gene_list = [str(introns_list[0])+str(exon_list[0])+str(introns_list[1])]
            Gene.total_length= len(Gene.full_gene_list[0])
            Gene.preintron = str(introns_list[0])
            Gene.postintron = str(introns_list[1])
            Gene.exon= str(exon_list[0])           
            if Gene.header not in gene_dict: 
                gene_dict[Gene.header]=[Gene.preintron, Gene.exon, Gene.postintron, Gene.full_gene_list, Gene.total_length]
        # print("gene dictionary", gene_dict)
    return 

genie()


motif_dict={}

#Dictionary of IPUAC degenerate bases (as chars)

degen_bases_dict={  "u" : "[tu]",
                    "y" : "[tcu]",
                    "n" : "[agtc]",
                    "U" : "[tu]",
                    "Y" : "[tcu]",
                    "N" : "[agtc]",
                    "a" : "[a]",
                    "c" : "[c]",
                    "g" : "[g]",
                    "t" : "[tu]",
                    } 

class Motifs:
    def __init__(self, motif, degen_introns, degen_exons, position):
            '''The Gene class'''

            self.motif = motif
            self.degen_introns = degen_introns
            self.degen_exons = degen_exons
            self.position = position

def motif_maker():
    with open(motif , "r") as mfh:
        for line in mfh:
            line = line.strip().lower()
            standard = line
            rmotif = ""
            for letter in line:
                rmotif += degen_bases_dict[letter]
            if standard not in motif_dict:
                motif_dict[standard]=rmotif 
        # print(motif_dict)
    return
    
motif_maker()

motif_position_dict={}

def motif_searcher():
        for l in range(0,4):
            v =  list(gene_dict.values())
            stringer = str(str(v[l][3][0]))
            stringer = stringer.lower()
            for j in range(0,4): 
                v = list(motif_dict.values())
                k = list(motif_dict.keys())
                query = v[j]
                mot = k[j]
                for m in re.finditer(query, stringer):
                    start = m.start()
                    mot = mot 
                    
                    if stringer not in motif_position_dict:
                        motif_position_dict[stringer]=[(start, mot)]
                    else:
                        motif_position_dict[stringer]+=[(start, mot)]
        # print(motif_position_dict)

motif_searcher()

## generating motif objects 

def getList(dict):
    return list(dict.keys())

header_list = getList(gene_dict)

motif_object_dict = {}

i =0
for key in motif_dict: 
    motif_dict[key]=color_pal_list[i]
    i+=1
    # print(motif_dict)

i = 0 
offset = 0

for i in range(0, len(header_list)):
    
    preintron = gene_dict.get(header_list[i])[0]
    exon = gene_dict.get(header_list[i])[1]
    postintron = gene_dict.get(header_list[i])[2] 
    full_gene = gene_dict.get(header_list[i])[3]
    full_gene = (full_gene[0])
    full_gene_lower = full_gene.lower()
    gene_length = gene_dict.get(header_list[i])[4]
    if i not in motif_object_dict:
        motif_object_dict[i]=[header_list[i], preintron, exon, postintron, gene_length, offset, motif_position_dict[full_gene_lower] ]  # motif_position_dict[i]
    else:
        motif_object_dict[i]+=[header_list[i], preintron, exon, postintron, gene_length, offset, motif_position_dict[full_gene_lower]]
        i = i + 1
    offset+= 250


i=0 
max_length = 0  

for i in range(0, len(motif_object_dict)):
    
    header = motif_object_dict[i][0]
    prei = motif_object_dict[i][1]
    ex = motif_object_dict[i][2]
    post = motif_object_dict[i][3]
    length=motif_object_dict[i][4]
    offset = motif_object_dict[i][5]
    motifs = motif_object_dict[i][6]
    # print(length) # works 
    j = 0 
    for j in range(0, len(motif_object_dict)):
        if (length > max_length) == True:
            max_length = length
            j+=1 
            # print(i, max_length) # works 
    i+=1
    

#cairo 

# set width and height based on the length of the longest gene as well as the number of genes to be viewed. 

width = int(1.25 * max_length)
height = len(motif_object_dict) * 250 #? 
num_objects = len(motif_object_dict)


surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, width, height)
context = cairo.Context(surface)
context.rectangle(0, 0, width, height) # making a surface-sized rectangle 
context.set_source_rgb(1, 1, 1) # making it white 
context.fill() # filling it in 

for i in range(0, len(motif_object_dict)): 

    header = motif_object_dict[i][0]
    prei = motif_object_dict[i][1]
    len_prei = len(motif_object_dict[i][1])
    ex = motif_object_dict[i][2]
    len_ex = len(motif_object_dict[i][2])
    post = motif_object_dict[i][3]
    len_post = len(motif_object_dict[i][3])
    length=motif_object_dict[i][4]
    offset = motif_object_dict[i][5]
    motifs = motif_object_dict[i][6]

    start_spot_x = (.1*width )
    start_spot_y = (height/num_objects)/2 + offset

    context.set_source_rgb(0,0,0) # change pen to black
    context.set_line_width(1) # fine point pen 

    context.move_to(start_spot_x,start_spot_y)        #(x,y) origin (0,0) is to .   .p left

    context.line_to(start_spot_x+length,start_spot_y)
    context.stroke()

    context.move_to(start_spot_x-10, start_spot_y-15)
    context.select_font_face("Lato", cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
    context.set_font_size(16)
    context.set_source_rgb(0,0,0)
    context.show_text(header)

    context.rectangle((start_spot_x+len_prei),(start_spot_y - ((height*0.0125)/2) ),(length-(len_prei + len_post)),height*0.0125)     
    #syntax for .rectangle is (beginning x, beginning y, width, height) so you just need the start and the dims
    
    context.fill()   

    #marking the motifs 

    #for loop to access each motif within each list
    i = 0 
    for i in range(0,len(motifs)):
        # print(motifs[i][0]) #works 
         
    #move to start_spot_x, start_spot_y, #moving to the origin of the gene
        
        motif_start_pos = motifs[i][0]
        spec_motif = motifs[i][1]
        motif_length = len(spec_motif)
        x = start_spot_x+motif_start_pos
        y = start_spot_y
 
        context.move_to(x,y)
 
        if spec_motif in motif_dict:
            color = "%.4f, %.4f, %.4f" % (motif_dict[spec_motif])
 
        color1=float(color[0:5])
        color2=float(color[8:13])
        color3=float(color[16:20])
 
        context.rectangle(x, y - ((height*0.0125)/2), motif_length, height*0.0125 )
        context.set_source_rgb(color1,color2,color3)
        context.fill()
 
  
        i+=1
            
    # making a ledgend 
    
    context.move_to(start_spot_x-10, start_spot_y+25)

    font = context.select_font_face("Lato", cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
    fontsize = context.set_font_size(16)
    font_color = context.set_source_rgb(0,0,0)

    context.select_font_face("Lato", cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
    context.set_font_size(16)
    context.set_source_rgb(0,0,0)
    context.show_text("Key:")

    x_offset = 50
    y_offset=20 
    i=0
    for spec_motif in motif_dict:
        move = context.move_to(start_spot_x+x_offset, start_spot_y+25)
        move
        x_offset+=40
        move 
        if spec_motif in motif_dict:
            color = "%.4f, %.4f, %.4f" % (motif_dict[spec_motif])
 
            color1=float(color[0:5])
            color2=float(color[8:13])
            color3=float(color[16:20])
 
        context.set_source_rgb(color1,color2,color3)
        context.show_text(f"{spec_motif}")
        # context.rectangle(start_spot_x+x_offset, start_spot_y+50- ((height*0.0125)/2), motif_length, height*0.0125)
        
        x_offset+=20
       
    context.move_to(start_spot_x ,start_spot_y ) # for the next motif object 
    prefix = args.filename.split(".")[0]
surface.write_to_png(f"{prefix}"+".png") 
surface.finish()