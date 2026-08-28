#!/usr/bin/env python3
import cgi
import cgitb; cgitb.enable()  # for troubleshooting
from os.path import abspath,dirname,isfile
import os
import subprocess
from string import Template
import sys
import re
from datetime import datetime

min_seq_len=3       # minimum sequence length
max_seq_len=10000   # maximum sequence length

outdir =abspath(dirname(__file__)+"/../output")
cif2pdb=abspath(dirname(__file__)+"/module/cif2pdb")
pdb2fasta=abspath(dirname(__file__)+"/module/pdb2fasta")
TMalign=abspath(dirname(__file__)+"/module/USalign")
bywhom=abspath(dirname(__file__)+"/../log/bywhom.txt")

# Added flex_args for advanced flexible alignment options
cmd_template=Template("cd $outdir;/usr/bin/timeout 1m $TMalign -mm ${mm} ${jobID}A.pdb ${jobID}B.pdb -byresi ${byresi} -mol ${mol} -atom \"${atom}\" ${flex_args} -o ${jobID}super")

html_template=Template('''<html>
<head><title>US-align $jobID</title></head>
<script type="text/javascript" src="/jmol/JSmol.min.js"></script>
<script type="text/javascript" src="/jmol/Jmol2.js"></script>

<body bgcolor=#f6f6f6>
<a href=..>[Back to server]</a>
<h3>Alignment result for job $jobID</h3>

<ul>
<pre>
$TMalign_output 
</pre>
</ul>

<h3>Superposition (Structure_1 in color and Structure_2 in white/black)</h3>
<ul>
<table>
<tr><td>
<script type="text/javascript">
jmolInitialize("/jmol/");
jmolSetAppletColor("#000000");
/* Injected jmol_color_cmds parsed from PyMOL script for dynamic hinge coloring */
jmolApplet(500, "load ${jobID}.pdb; calculate structure rama; select all; cartoons; spacefill off; wireframe off; backbone off; select :B and not ligand; cartoon off; trace on; ${jmol_color_cmds} select :A and ligand; color cyan; select :B and ligand; color pink; set spin x 10; set spin y 10; spin off");
</script>
</td></tr>
<tr><td>
<script>
jmolButton("Reset","Reset to initial orientation");
jmolCheckbox("spin on", "spin off", "Spin On/Off", false);
jmolCheckbox("set antialiasDisplay true",
    "set antialiasDisplay false", "High quality", false);
</script>
</td></tr>
<tr><td>
<script>
jmolButton("write PNGJ USalign.png","Save image");
jmolCheckbox("color background white; select :B and not ligand; color black",
    "color background black; select :B and not ligand; color white", "White background", false);
jmolCheckbox("select :B and not ligand; cartoon off; trace on",
    "select :B and not ligand; cartoons; trace off",
    "Show Structure_2 in trace", true);
</script>
</td></tr>
<tr><td>
<script>
jmolCheckbox("select ligand; wireframe 0.25; spacefill 0.5; select solvent; spacefill 0.25;",
    "select all; wireframe off; spacefill off;", "Show ligand", false);
</script>
(water of Structure_1 and Structure_2 in blue and red, respectively;<br>
non-water ligand of Structure_1 and Structure_2 in cyan and pink, respectively).
</td></tr>
</table>
</ul>

<h3>Result download</h3>
<ul>
<li>Structure_1 (before superposition): <a href=${jobID}A.pdb target=_blank>
${jobID}A.pdb</a>.</li>
<li>Structure_2 (before superposition): <a href=${jobID}B.pdb target=_blank>
${jobID}B.pdb</a>.</li>
<li>Structure_1 superposed to Structure_2 (with Structure_2): 
<a href=${jobID}.pdb target=_blank>
${jobID}.pdb</a>.</li>
<li>Structure_1 superposed to Structure_2 (without Structure_2): 
<a href=${jobID}super.pdb target=_blank>
${jobID}super.pdb</a>.</li>
All the above files will be removed from the server after 2 days.
</ul>

<a href=..>[Back to server]</a>
<p>
</body>
<html/>
''')

def parse_pml_to_jmol(pml_filename):
    ''' Parse the generated _all_atm.pml PyMOL script to generate Jmol color commands '''
    if not isfile(pml_filename):
        # Fallback to default coloring if PML is not found
        return "select :A; color blue; select :B; color red;"
    
    jmol_cmds = []
    selections = {}
    
    # Map common PyMOL gray scales to Jmol RGB arrays
    color_map = {
        "gray70": "[178,178,178]",
        "gray50": "[127,127,127]",
        "gray30": "[76,76,76]",
        "grey"  : "[128,128,128]"
    }
    
    try:
        fp = open(pml_filename, 'r')
        for line in fp:
            line = line.strip()
            
            # Extract selections (e.g. select hinge_0_str1, structure1 and c. A and ( i. 20 or i. 21 ...))
            if line.startswith("select"):
                parts = line.split(',', 1)
                if len(parts) > 1:
                    sel_name = parts[0].split()[1].strip()
                    condition = parts[1]
                    
                    chain = ":A" if "structure1" in condition else ":B"
                    
                    # Extract all residue indices using regex
                    indices = re.findall(r'i\.\s*(-?\d+)', condition)
                    if indices:
                        selections[sel_name] = (chain, indices)
                        
            # Extract colors (e.g. color gray70, structure1 or color blue, hinge_0_str1)
            elif line.startswith("color"):
                parts = line.split(',', 1)
                if len(parts) > 1:
                    c_name = parts[0].split()[1].strip()
                    target = parts[1].strip()
                    
                    jmol_color = color_map.get(c_name, c_name)
                    
                    if target == "structure1":
                        jmol_cmds.append("select :A; color %s;" % jmol_color)
                    elif target == "structure2":
                        jmol_cmds.append("select :B; color %s;" % jmol_color)
                    elif target in selections:
                        chain, indices = selections[target]
                        idx_str = ",".join(indices)
                        jmol_cmds.append("select %s and (%s); color %s;" % (chain, idx_str, jmol_color))
        fp.close()
    except:
        print("ERROR! cannot parse "+pml_filename)
        
    if not jmol_cmds:
        return "select :A; color blue; select :B; color red;"
        
    return " ".join(jmol_cmds)

def format_alignment_output(output, pml_filename):
    '''
    Parse the US-align output sequence block, mapping colors from the PML file to 
    Structure 1, hinge-matching line, and Structure 2 without bolding to avoid misalignment.
    '''
    # Default PyMOL colors typically used by US-align hinges
    colors = {
        '0': 'blue', '1': 'red', '2': 'green', '3': 'yellow',
        '4': 'magenta', '5': 'cyan', '6': 'orange', '7': 'springgreen',
        '8': 'teal', '9': 'purple', 'a': 'pink', 'b': 'brown', 'c': 'gray',
        ':': 'blue', '.': 'blue'
    }
    # Convert specific PyMOL colors to CSS-compatible hex values
    css_map = {
        "gray70": "#b2b2b2", "gray50": "#7f7f7f", 
        "gray30": "#4c4c4c", "grey": "grey",
        "spring": "springgreen"
    }
    
    # Dynamically read exact colors from the .pml file if available
    try:
        if isfile(pml_filename):
            with open(pml_filename, 'r') as fp:
                for line in fp:
                    line = line.strip()
                    if line.startswith('color '):
                        parts = line.split(',', 1)
                        if len(parts) == 2:
                            cname = parts[0].split()[1].strip()
                            target = parts[1].strip()
                            m = re.search(r'hinge_([0-9a-zA-Z])', target)
                            if m:
                                colors[m.group(1)] = cname
    except:
        pass

    # Replace carriage returns for safer manipulation, then split by newlines
    lines = output.replace('\r', '').split('\n')
    start_idx = -1
    
    # Locate the starting line of the sequence alignment block
    for i, line in enumerate(lines):
        if "denote different aligned fragment pairs" in line or \
           "denotes residue pairs" in line or \
           "denotes aligned residue pairs" in line:
            start_idx = i + 1
            break

    if start_idx != -1:
        align_indices = []
        for i in range(start_idx, len(lines)):
            # Skip empty lines or the Total CPU time comment line at the end
            if lines[i].strip() and not lines[i].startswith("#"):
                align_indices.append(i)
                
        # Process every block of 3 lines (Seq 1, Align match/hinge line, Seq 2)
        for b in range(0, len(align_indices) - 2, 3):
            idx1 = align_indices[b]
            idx_mid = align_indices[b+1]
            idx2 = align_indices[b+2]
            
            seq1 = lines[idx1]
            mid = lines[idx_mid]
            seq2 = lines[idx2]
            
            # Ensure equal string lengths for smooth iteration over sequences
            max_len = max(len(seq1), len(mid), len(seq2))
            seq1 = seq1.ljust(max_len)
            mid = mid.ljust(max_len)
            seq2 = seq2.ljust(max_len)
            
            n_seq1, n_mid, n_seq2 = "", "", ""
            for c1, cm, c2 in zip(seq1, mid, seq2):
                # A character other than space or '-' in the middle line indicates an aligned residue
                if cm != ' ' and cm != '-':
                    col = css_map.get(colors.get(cm, 'blue'), colors.get(cm, 'blue'))
                    n_seq1 += '<span style="color:%s">%s</span>' % (col, c1)
                    n_mid  += '<span style="color:%s">%s</span>' % (col, cm)
                    n_seq2 += '<span style="color:%s">%s</span>' % (col, c2)
                else:
                    n_seq1 += c1
                    n_mid  += cm
                    n_seq2 += c2
                    
            lines[idx1] = n_seq1
            lines[idx_mid] = n_mid
            lines[idx2] = n_seq2
            
    return '\n'.join(lines)


def read_structure_as_chain(filename,chainID):
    ''' read structure '''
    txt=''
    fp=open(filename,'r')
    for line in fp.read().splitlines():
        if line.startswith("ATOM  ") or line.startswith("HETATM"):
            txt+=line[:21]+chainID+line[22:]+'\n'
    fp.close()
    return txt

def make_display_structure(outdir,jobID):
    ''' combine $outdir/${jobID}super.pdb and $outdir/${jobID}B.pdb into
    $outdir/${jobID}.pdb'''
    fp=open(outdir+'/'+jobID+".pdb",'w')
    if isfile(outdir+'/'+jobID+"super.pdb"):
        fp.write(read_structure_as_chain(outdir+'/'+jobID+"super.pdb",'A'))
    else:
        fp.write(read_structure_as_chain(outdir+'/'+jobID+"A.pdb",'A'))
    fp.write(read_structure_as_chain(outdir+'/'+jobID+"B.pdb",'B'))
    fp.close()
    return

def check_structure(pdbin,ter='3',atom="auto",mol="auto"):
    '''check if the user input structure is correct
    if it is correct, return the sequence
    if it is incorrect, print the error message and die
    '''
    cmd = "%s - -ter %s -atom \"%s\" -mol %s"%(pdb2fasta,ter,atom,mol)
    p = subprocess.Popen(cmd,
        stdout=subprocess.PIPE, stdin=subprocess.PIPE, shell=True)
    fastatxt,stderr = p.communicate(input=pdbin.encode())
    fastatxt=fastatxt.decode()
    sequence=''
    for line in fastatxt.splitlines():
        sequence+=line
    if len(sequence)<=min_seq_len:
        print("ERROR! Too few residues: %d residues (<=%d)"%(
            len(sequence),min_seq_len))
        exit()
    if len(sequence)>max_seq_len:
        print("ERROR! Too many residues: %d residues (>%d). Please download the standalone program for large input."%(
            len(sequence),max_seq_len))
        exit()

    if "\nloop_" in pdbin: # PDBx/mmCIF format input
        p = subprocess.Popen(cif2pdb+" -",
            stdout=subprocess.PIPE, stdin=subprocess.PIPE, shell=True)
        pdbin,stderr = p.communicate(input=pdbin.encode())
        pdbin=pdbin.decode()
        print(stderr)

    pdbtxt=''
    for line in pdbin.splitlines():
        if len(sequence)>0:
            if int(ter)>=1 and line.startswith("END"):
                break # remove multiple model
        if not line.startswith("ATOM  ") and \
           not line.startswith("HETATM") and \
           not line.startswith("CONECT") and \
           not line.startswith("MODEL") and \
           not line.startswith("END") and \
           not line.startswith("TER"):
            continue
        pdbtxt+=line+'\n'
    return sequence,pdbtxt

def get_cgi_field():
    '''get cgi input'''
    form = cgi.FieldStorage()
    pdbA=form.getfirst("pdb71_txt",'')
    pdbB=form.getfirst("pdb72_txt",'')
    if not pdbA.strip():
        pdbA=form.getfirst("pdb71_file",'').strip()
        if len(pdbA):
            pdbA=pdbA.decode()
    if not pdbB.strip():
        pdbB=form.getfirst("pdb72_file",'').strip()
        if len(pdbB):
            pdbB=pdbB.decode()
    byresi=form.getfirst("byresi","0")
    mol=form.getfirst("mol","auto")
    atom=form.getfirst("atom","auto")
    mm=form.getfirst("mm","7")
    
    # Advanced flexible alignment arguments
    fast=form.getfirst("fast","0")
    hinge=form.getfirst("hinge","5")
    afp=form.getfirst("afp","0")
    TMpass=form.getfirst("TMpass","0.85")
    
    return str(pdbA),str(pdbB),byresi,mol,atom,mm,hinge,afp,TMpass,fast

def write_bywhom(jobID,sequenceA,sequenceB):
    ip=str(os.getenv("REMOTE_ADDR"))
    fp=open(bywhom,'a')
    fp.write('\t'.join([jobID,ip,str(len(sequenceA)),str(len(sequenceB)),
        str(datetime.now()),sequenceA,sequenceB])+'\n')
    fp.close()
    return

if __name__=="__main__":
    print("Content-type: text/html\n")
    print("<head><title>US-align</title></head>\n")
    pdbA,pdbB,byresi,mol,atom,mm,hinge,afp,TMpass,fast=get_cgi_field()

    if len(sys.argv)>1:
        # CGI script does not normally receive command line argument.
        # The command line arguments are for debug purpose only
        for arg in sys.argv[1:]:
            if arg.startswith("-pdb1="):
                fp=open(arg[len("-pdb1="):],'r')
                pdbA=fp.read()
                fp.close()
            elif arg.startswith("-pdb2="):
                fp=open(arg[len("-pdb2="):],'r')
                pdbB=fp.read()
                fp.close()
            elif arg.startswith("-byresi="):
                byresi=arg[len("-byresi="):]
            elif arg.startswith("-mol="):
                mol=arg[len("-mol="):]
            elif arg.startswith("-atom="):
                atom=arg[len("-atom="):]
            elif arg.startswith("-mm="):
                mm=arg[len("-mm="):]
            elif arg.startswith("-fast="):
                fast=arg[len("-fast="):]
            elif arg == "-fast":
                fast="1"
            elif arg.startswith("-hinge="):
                hinge=arg[len("-hinge="):]
            elif arg.startswith("-afp="):
                afp=arg[len("-afp="):]
            elif arg.startswith("-TMpass="):
                TMpass=arg[len("-TMpass="):]
            else:
                print("ERROR! No such option "+arg)

    sequenceA,pdbA=check_structure(pdbA,atom=atom,mol=mol)
    sequenceB,pdbB=check_structure(pdbB,atom=atom,mol=mol)

    jobID=subprocess.Popen("date +%N",shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ).communicate()[0].strip()
    jobID="US"+jobID.decode()
    
    fp=open(outdir+'/'+jobID+"A.pdb",'w')
    fp.write(pdbA)
    fp.close()
    fp=open(outdir+'/'+jobID+"B.pdb",'w')
    fp.write(pdbB)
    fp.close()
        
    print('<meta http-equiv="refresh" content="0.1;url=../output/%s.html">'%jobID)

    # Process flexible alignment specific flags
    flex_args = ""
    if fast == "1":
        flex_args += "-fast "
    if hinge and hinge.strip().isdigit():
        flex_args += "-hinge " + hinge.strip() + " "
    if afp == "1":
        flex_args += "-afp "
        if TMpass:
            flex_args += "-TMpass " + TMpass.strip() + " "

    cmd=cmd_template.substitute(outdir=outdir,TMalign=TMalign,jobID=jobID,
        byresi=byresi,mol=mol,atom=atom,mm=mm,flex_args=flex_args)
    TMalign_output,stderr=subprocess.Popen(cmd,shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    TMalign_output=TMalign_output.decode()
    stderr=stderr.decode()

    # Parse PyMOL script for dynamic Jmol coloring
    pml_filename = outdir + '/' + jobID + "super_all_atm.pml"
    if not isfile(pml_filename):
        pml_filename = outdir + '/' + jobID + "super.pml" # fallback
        
    jmol_color_cmds = parse_pml_to_jmol(pml_filename)

    # ==========================================
    # Call the newly added formatting function
    # to apply colors to the alignment sequences
    # ==========================================
    TMalign_output = format_alignment_output(TMalign_output, pml_filename)

    make_display_structure(outdir, jobID)
    fp=open(outdir+'/'+jobID+".html",'w')
    fp.write(html_template.substitute(
        jobID=jobID, 
        TMalign_output=TMalign_output.strip(),
        jmol_color_cmds=jmol_color_cmds
    ))
    fp.write(stderr)
    fp.close()

    # failure to write bywhom will not affect user output display
    write_bywhom(jobID,sequenceA,sequenceB)
