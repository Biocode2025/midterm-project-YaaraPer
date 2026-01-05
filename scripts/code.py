def DNA_RNA_Cod(seq):
  '''
  הפונקציה דואגת שהאותיות תהיינה אחידות (אותיות גדולות) והופכת את רצף ה- DNA המקודד לרצף RNA.
  מקבלת: seq.
  מחזירה: RNA_seq.
  '''
  RNA_seq=""

  seq = seq.upper()
  
  rna_line = seq.replace("T","U")
  RNA_seq += rna_line
    
  return RNA_seq
#------------------------------------------------

def Read_dict(fl):
  '''
  הפונקציה קוראת לתוךdictionary  את המיפוי בין הקודונים לחומצות אמינו מהקובץ.
  מקבלת: fl.
  מחזירה: RNA_codon_table.
  '''
  global RNA_codon_table

  for line in fl:
    line = line.rstrip('\r\n')
    line_list = line.split()
    
    codon_key = line_list[0]
    amino_acids_value = line_list[1]
    
    RNA_codon_table[codon_key] = amino_acids_value
    
  return RNA_codon_table
#------------------------------------------------

def RNA_prot(seq):
  '''
  הפונקציה מתרגמת את רצף ה- RNA לרצף חלבון.
  מקבלת: seq.
  מחזירה: protein_seq.
  '''
  protein_seq = ""
  
  for i in range(0, len(seq), 3):
    codon = seq[i:i+3]
    if len(codon) < 3:
      break
    
    curr_amino_acid = RNA_codon_table.get(codon,"-")
    protein_seq += curr_amino_acid
  
  return protein_seq
#------------------------------------------------

def get_operon_seq(start, stop, full_genom):
  '''
  הפונקציה מחלצת מהגנום כולו אופרונים על פי מיקום התחלה וסוף שהיא מקבלת.
  מקבלת: start, stop, full_genom
  מחזירה: operon
  '''
  operon = full_genom[start:stop + 1]
  return operon
#------------------------------------------------

def get_sub_gene(operon_seq):
  '''
  הפעולה מחלצת מהאופרון גנים. במידה והרצף שהתקבל (כל מה שבא אחרי הגן הנוכחי) מספיק ארוך, ניתן לקרוא אותו לפי קודונים בצורה תקינה- הפונקציה תחזיר את הגן ושאר הרצף אחרי הגן.
  אם הרצף שהתקבל לא מספיק ארוך הפונקציה תחזיר מחרוזת ריקה.
  מקבלת: operon_seq
  מחזירה: operon_seq[starting_point:ending_point], operon_seq[ending_point:] או ""
  '''

  starting_point = 0
  ending_point = 0
  operon_len = len(operon_seq)
  
  if len(operon_seq) >= 6:
    for i in range(0, operon_len, 3):
      codon = operon_seq[i:i+3]
      
      if codon == "ATG":
        starting_point = i
        
        for h in range(starting_point + 3, operon_len, 3):
          codon = operon_seq[h:h+3]
          
          if codon == "TAA" or codon == "TAG"  or codon == "TGA":
            ending_point = h + 3
            
            return operon_seq[starting_point:ending_point], operon_seq[ending_point:]
  
  return ""
#------------------------------------------------
def split_line(ln):
  '''
  הפונקציה מפצלת שורה בטבלה לכמה משתנים שונים ורשימה אחת.
  מקבלת: ln
  מחזירה: int(starting_point), int(eding_point), strand, int(number_of_genes), ln_list
  '''
  
  starting_point, eding_point, strand, number_of_genes, temp_str = ln.split("\t", 4)
    
  cnt = 0
  new_ln = ""
  
  ln_list = temp_str.split(",")
  return int(starting_point), int(eding_point), strand, int(number_of_genes), ln_list
  
  

# תוכנית ראשית.
global RNA_codon_table
RNA_codon_table = {}

# הגדרת משתנים
ecoli_genom = ""
number_of_operon = 0
Num_genes = 0
cnt = 0
all_genes = 0
percentages = 0

# פתיחת קבצים לקריאה
ecoli_genom_file = open('data/ecoli.full.fasta.txt', 'r')
codon_file = open('data/codon_AA.txt', 'r')
Read_dict(codon_file)
operons_table_file = open('data/e.coli_operons.txt','r')

# פתיחת קובץ לכתיבה.
proteines_file = open('results/proteines.txt', 'w')


# הלולאה קוראת את רצף ה- DNA של החיידק מהקובץ.
for line in ecoli_genom_file:
  line = line.rstrip('\r\n')
  if line == "":
    continue
  
  # רצף ה DNA מופיע בשורות שאינן מתחילות בסימן "<" לכן "נדלג" על שורה זו
  if line[0] == ">":
    continue
  
  ecoli_genom += line
  

# לולאה חיצונית: קוראת את הקובת שבו יש את הטבלה ומשתמשת בנתונים שבטבלה.
for line in operons_table_file:
  line = line.rstrip('\n')
  
  # התוכנית "מדלגת על השורה הראשונה של הטבלה."
  if line[0] != 'S':
    # בכל איטרציה של הלולאה משתנה האופרון לכן בתחילת הלולאה "נמחק" את רצף האופרון הקודם.
    curr_operon = ""

    # קריאה לפונקציה שמחלקת שורה בטבלה למספר משתנים.
    Start, Stop, Strand, Num_genes, Genes_names = split_line(line)
    if Strand == "+":
      number_of_operon = number_of_operon + 1
      
      # קריאה לפונקציה שמחזירה את רצף האופרון הנוכחי, לפי הנתונים בטבלה.
      curr_operon = get_operon_seq(Start, Stop, ecoli_genom)
      
      # כתיבה בקובץ.
      proteines_file.write("Operon number ")
      proteines_file.write(str(number_of_operon))
      proteines_file.write(":"  + '\n')
      
      # לולאה פנימית: בלולאה זו מוצאים את הגנים באופרון, "מתרגמים" אותם וכותבים את התוצאות הסופיות בקובץ.
      for i in range(Num_genes):
        
        # הוראת תנאי, היא מונעת עבודה על מחרוזת ריקה.
        if (get_sub_gene(curr_operon) != ""):
          # קריאה לפונקציות, כדי "לחלץ" גן מהאופרון ולתרגם גן זה.
          gene, curr_operon = get_sub_gene(curr_operon)
          gene_as_RNA = DNA_RNA_Cod(gene)
          protein = RNA_prot(gene_as_RNA)
          
          # חיבור האורכים של כל הגנים.
          all_genes = all_genes + len(gene)
          
          # כתיבה בקובץ.
          proteines_file.write(Genes_names[i].strip(" "))
          proteines_file.write(": ")
          proteines_file.write(protein + '\n')
          
          gene = ""
        
        else:
          break
      proteines_file.write('\n')
      cnt = cnt + 1  
    

# חישוב כמה אחוזים מהגנום של החיידק מקודדים לחלבונים.
percentages = (all_genes / len(ecoli_genom)) * 100
print("Percentages of the bacterial genome that encode proteins: %.2f%%" %percentages)
 
# סגירת הקבצים
proteines_file.close()
ecoli_genom_file.close()
codon_file.close()
operons_table_file.close()