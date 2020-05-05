# annotate_v
Antibody Annotation - A simple to use Python script annotating VH or VL sequences of an antibody using **Kabat**, **Chothia**, **Contact**, or **IMGT** schemes. It utilizes the REST API interface of [Abnum](http://www.bioinf.org.uk/abs/abnum/) from Dr Andrew Martin's group at UCL

# Description

User provides a single-letter amino acid sequence of the VH or VL chain of an antibody, and specifies an annotation scheme (Kabat, Chothia, or Contact, or IMGT). The script sends a request to [Abnum](http://www.bioinf.org.uk/abs/abnum/), which returns a string with each residue matched to its number. The script then identifies FR and CDR regions using definitions outlined [here](http://www.bioinf.org.uk/abs/info.html#kabatnum)<sup>1</sup>. It prints out the amino acid sequence of the FRs and CDRs, and returns a `list` of 2 `dict`. The first `dict` consists of `region: seq` pairs. The second `dict` consists of `number:residue` pairs. 

<sup>1</sup>Definitions for Kabat, Chothia, Contact, and IIMGT are same as listed in the table, except that for IMGT, H-CDR2 is defined as **H51-H57** in this script, as opposed to of **H51-H56** in the table. This slight revision generates result that matches that from [IMGT website](http://www.imgt.org/)

# Dependencies
- Imports `requests` module

# Simple to use

## First assign variables and create an instance of *class* `annotate`

`aaseq`: *STRING* amino acid sequence, **single-letter** coded, needs to be **complete VH OR VL** sequence. Upper/lower case. 

`scheme`:*STRING* annotation scheme, can be one of the following: **"kabat", "chothia", "contact", "imgt"**. Must be **lowercase**

```python
aaseq="EIVLTQSPAIMSASPGERVTMTCSASSGVNYMHWYQQKPGTSPRRWIYDTSKLASGVPARFSGSGSGTDYSLTISSMEPEDAATYYCHQRGSYTFGGGTKLEIK"
scheme="chothia"
seq=annotate(aaseq,scheme)
```

## Then invoke the `retrieve()` method

```python
result=seq.retrieve()
print(result[0]) #prints the first dict (region vs seq)
print(result[1]) #prints the second dict (number vs residue)

```

## Example output

```
Annotation in Chothia scheme:
L-FR1:   EIVLTQSPAIMSASPGERVTMTC
L-CDR1:  SASSGVNYMH
L-FR2:   WYQQKPGTSPRRWIY
L-CDR2:  DTSKLAS
L-FR3:   GVPARFSGSGSGTDYSLTISSMEPEDAATYYC
L-CDR3:  HQRGSYT
L-FR4:   FGGGTKLEIK
{'L-FR1': 'EIVLTQSPAIMSASPGERVTMTC', 'L-CDR1': 'SASSGVNYMH', 'L-FR2': 'WYQQKPGTSPRRWIY', 'L-CDR2': 'DTSKLAS', 'L-FR3': 'GVPARFSGSGSGTDYSLTISSMEPEDAATYYC', 'L-CDR3': 'HQRGSYT', 'L-FR4': 'FGGGTKLEIK'}
{'L1': 'E', 'L2': 'I', 'L3': 'V', 'L4': 'L', 'L5': 'T', 'L6': 'Q', 'L7': 'S', 'L8': 'P', 'L9': 'A', 'L10': 'I', 'L11': 'M', 'L12': 'S', 'L13': 'A', 'L14': 'S', 'L15': 'P', 'L16': 'G', 'L17': 'E', 'L18': 'R', 'L19': 'V', 'L20': 'T', 'L21': 'M', 'L22': 'T', 'L23': 'C', 'L24': 'S', 'L25': 'A', 'L26': 'S', 'L27': 'S', 'L28': 'G', 'L29': 'V', 'L30': 'N', 'L32': 'Y', 'L33': 'M', 'L34': 'H', 'L35': 'W', 'L36': 'Y', 'L37': 'Q', 'L38': 'Q', 'L39': 'K', 'L40': 'P', 'L41': 'G', 'L42': 'T', 'L43': 'S', 'L44': 'P', 'L45': 'R', 'L46': 'R', 'L47': 'W', 'L48': 'I', 'L49': 'Y', 'L50': 'D', 'L51': 'T', 'L52': 'S', 'L53': 'K', 'L54': 'L', 'L55': 'A', 'L56': 'S', 'L57': 'G', 'L58': 'V', 'L59': 'P', 'L60': 'A', 'L61': 'R', 'L62': 'F', 'L63': 'S', 'L64': 'G', 'L65': 'S', 'L66': 'G', 'L67': 'S', 'L68': 'G', 'L69': 'T', 'L70': 'D', 'L71': 'Y', 'L72': 'S', 'L73': 'L', 'L74': 'T', 'L75': 'I', 'L76': 'S', 'L77': 'S', 'L78': 'M', 'L79': 'E', 'L80': 'P', 'L81': 'E', 'L82': 'D', 'L83': 'A', 'L84': 'A', 'L85': 'T', 'L86': 'Y', 'L87': 'Y', 'L88': 'C', 'L89': 'H', 'L90': 'Q', 'L91': 'R', 'L92': 'G', 'L93': 'S', 'L96': 'Y', 'L97': 'T', 'L98': 'F', 'L99': 'G', 'L100': 'G', 'L101': 'G', 'L102': 'T', 'L103': 'K', 'L104': 'L', 'L105': 'E', 'L106': 'I', 'L107': 'K'}
```

# Limitations
- Relies on internet connection to [Abnum](http://www.bioinf.org.uk/abs/abnum/) website
- Currently it can only annotate using Kabat, Chothia, or Contact, or IMGT scheme. One scheme each time.
- Incomplete VH or VL sequence might not be annotated 
