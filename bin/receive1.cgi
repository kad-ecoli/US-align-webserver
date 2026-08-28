#!/usr/bin/env python3
import cgi
import cgitb; cgitb.enable()  # for troubleshooting
from os.path import abspath,dirname,isfile
import os
import subprocess
from string import Template
import sys
from datetime import datetime

min_seq_len=3       # minimum sequence length
max_seq_len=10000   # maximum sequence length
max_chain_num=62    # maximum number of chains

outdir =abspath(dirname(__file__)+"/../output")
cif2pdb=abspath(dirname(__file__)+"/module/cif2pdb")
pdb2fasta=abspath(dirname(__file__)+"/module/pdb2fasta")
TMalign=abspath(dirname(__file__)+"/module/USalign")
bywhom=abspath(dirname(__file__)+"/../log/bywhom.txt")

cmd_template=Template("cd $outdir;/usr/bin/timeout 1m $TMalign -mm 1 ${jobID}A.pdb ${jobID}B.pdb -ter ${ter} -mol ${mol} -atom \"${atom}\" -o ${jobID}super")

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

<h3>Superposition (Structure_1 in blue and Structure_2 in red)</h3>
<ul>
<table>
<tr><td>
<script type="text/javascript">
jmolInitialize("/jmol/");
jmolSetAppletColor("#000000");
jmolApplet(500, "load ${jobID}.pdb; calculate structure rama; select all; cartoons; spacefill off; wireframe off; backbone off; select :A; color blue; select :B; color red; select :A and ligand; color cyan; select :B and ligand; color pink; set spin x 10; set spin y 10; spin off"); 
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
jmolCheckbox("color background white",
    "color background black", "White background", false);
jmolCheckbox("select :B and not ligand; cartoon off; trace on",
    "select :B and not ligand; cartoons; trace off",
    "Show Structure_2 in trace", false);
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

def check_structure(pdbin,ter='1',atom="auto",mol="auto"):
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
    header_list=[]
    for line in fastatxt.splitlines():
        if line.startswith('>'):
            header_list.append(line)
        else:
            sequence+=line
    if len(header_list)>max_chain_num:
        print("ERROR! Too many chains: %d chains (>%d). Please download the standalone program for large input."%(
            len(header_list),max_chain_num))
        exit()
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
    pdbA=form.getfirst("pdb11_txt",'')
    pdbB=form.getfirst("pdb12_txt",'')
    if not pdbA.strip():
        pdbA=form.getfirst("pdb11_file",'').strip()
        if len(pdbA):
            pdbA=pdbA.decode()
    if not pdbB.strip():
        pdbB=form.getfirst("pdb12_file",'').strip()
        if len(pdbB):
            pdbB=pdbB.decode()
    ter=form.getfirst("ter",'1')
    mol=form.getfirst("mol","auto")
    atom=form.getfirst("atom","auto")
    return str(pdbA),str(pdbB),ter,mol,atom

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
    pdbA,pdbB,ter,mol,atom=get_cgi_field()

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
            elif arg.startswith("-ter="):
                ter=arg[len("-ter="):]
            elif arg.startswith("-mol="):
                ter=arg[len("-mol="):]
            elif arg.startswith("-atom="):
                ter=arg[len("-atom="):]
            else:
                print("ERROR! No such option "+arg)

    sequenceA,pdbA=check_structure(pdbA,ter,atom=atom,mol=mol)
    sequenceB,pdbB=check_structure(pdbB,ter,atom=atom,mol=mol)

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

    cmd=cmd_template.substitute(outdir=outdir,TMalign=TMalign,jobID=jobID,
        ter=ter,mol=mol,atom=atom)
    TMalign_output,stderr=subprocess.Popen(cmd,shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    TMalign_output=TMalign_output.decode()
    stderr=stderr.decode()

    make_display_structure(outdir, jobID)
    fp=open(outdir+'/'+jobID+".html",'w')
    fp.write(html_template.substitute(
        jobID=jobID, TMalign_output=TMalign_output.strip()))
    fp.write(stderr)
    fp.close()


    # failure to write bywhom will not affect user output display
    write_bywhom(jobID,sequenceA,sequenceB)
