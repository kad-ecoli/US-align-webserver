#!/usr/bin/env python3
import cgi
import cgitb; cgitb.enable()  # for troubleshooting
from os.path import abspath,dirname,isfile
import os
import subprocess
from string import Template
import sys
from datetime import datetime
import zipfile

min_seq_len=3       # minimum sequence length
max_seq_len=10000   # maximum sequence length
max_chain_num=10    # maximum number of chains

outdir =abspath(dirname(__file__)+"/../output")
cif2pdb=abspath(dirname(__file__)+"/module/cif2pdb")
pdb2fasta=abspath(dirname(__file__)+"/module/pdb2fasta")
TMalign=abspath(dirname(__file__)+"/module/USalign")
bywhom=abspath(dirname(__file__)+"/../log/bywhom.txt")

chainID_list="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"

cmd_template=Template("cd $outdir;/usr/bin/timeout 1m $TMalign -mm 4 -dir ${outdir}/ ${jobID}.list -suffix .pdb -mol ${mol} -atom \"${atom}\" -o ${jobID}super")

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

<h3>Superposition (Different Structure in different color)</h3>
<ul>
<table>
<tr><td>
<script type="text/javascript">
jmolInitialize("/jmol/");
jmolSetAppletColor("#000000");
jmolApplet(500, "load ${jobID}.pdb; calculate structure rama; select all; cartoons; spacefill off; wireframe off; backbone off; color chain; set spin x 10; set spin y 10; spin off"); 
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
jmolCheckbox("select not ligand; cartoon off; trace on",
    "select not ligand; cartoons; trace off",
    "Show structures in trace", false);
</script>
</td></tr>
<tr><td>
<script>
jmolCheckbox("select ligand; wireframe 0.25; spacefill 0.5; select solvent; spacefill 0.25;",
    "select all; wireframe off; spacefill off;", "Show ligand", false);
</script>
</td></tr>
</table>
</ul>

<h3>Result download</h3>
<ul>
${download_txt}
<li>Structures superposed to each other: <a href=${jobID}.pdb target=_blank>${jobID}.pdb</a>. (Chains in the superposed structure are sequentially re-assigned chain IDs of A, B, C, ... , J)</li>
All the above files will be removed from the server after 2 days.
</ul>

<a href=..>[Back to server]</a>
<p>
</body>
<html/>
''')

download_template=Template('''
<li>Structure_${I} (before superposition): <a href=${jobID}${I}.pdb target=_blank>${jobID}${I}.pdb</a>.</li>
''')

def read_structure(filename):
    ''' read structure '''
    txt=''
    fp=open(filename,'r')
    for line in fp.read().splitlines():
        if line.startswith("ATOM  ") or line.startswith("HETATM"):
            txt+=line+'\n'
    fp.close()
    return txt

def make_display_structure(outdir,jobID):
    ''' combine $outdir/${jobID}super.pdb and $outdir/${jobID}B.pdb into
    $outdir/${jobID}.pdb'''
    fp=open(outdir+'/'+jobID+".pdb",'w')
    for i in range(max_chain_num):
        filename=outdir+'/'+jobID+"super.%d.pdb"%i
        if isfile(filename):
            fp.write(read_structure(filename))
    fp.close()
    return

def check_structure(pdbin,ter='1',atom="auto",mol="auto",chainID=''):
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
        if len(chainID) and (line.startswith('ATOM  ') or line.startswith('HETATM')):
            if len(line)>22:
                line=line[:21]+chainID[-1]+line[22:]
        pdbtxt+=line+'\n'
    return sequence,pdbtxt

def get_cgi_field():
    '''get cgi input'''
    form = cgi.FieldStorage()
    pdb_list = []
    for i in range(max_chain_num):
        pdbA=form.getfirst("pdb4%d_txt"%(i+1),'')
        if not pdbA.strip():
            pdbA=form.getfirst("pdb4%d_file"%(i+1),'').strip()
            if len(pdbA):
                pdbA=pdbA.decode()
        if pdbA.strip():
            pdb_list.append(str(pdbA))
    mol=form.getfirst("mol","auto")
    atom=form.getfirst("atom","auto")

    fileitem = form['pdbzip_file']
    # Check if the file was uploaded successfully
    if fileitem.filename:
        # Open the zip file in memory
        zip_ref=zipfile.ZipFile(fileitem.file, 'r')
        # Process the files in the zip file
        for file_info in zip_ref.infolist():
            if not file_info.filename.endswith(".pdb") and \
               not file_info.filename.endswith(".ent") and \
               not file_info.filename.endswith(".cif"):
                continue
            # Extract the file
            fp=zip_ref.open(file_info)
            pdb_list.append(fp.read().decode())
            fp.close()
            if len(pdb_list)>=10:
                break
        zip_ref.close()
    return pdb_list,mol,atom

def write_bywhom(jobID,seq_len,sequence):
    ip=str(os.getenv("REMOTE_ADDR"))
    fp=open(bywhom,'a')
    fp.write('\t'.join([jobID,ip,seq_len,str(datetime.now()),sequence])+'\n')
    fp.close()
    return

if __name__=="__main__":
    print("Content-type: text/html\n")
    print("<head><title>US-align</title></head>\n")
    pdb_list,mol,atom=get_cgi_field()

    if len(sys.argv)>1:
        # CGI script does not normally receive command line argument.
        # The command line arguments are for debug purpose only
        for arg in sys.argv[1:]:
            if arg.startswith("-pdb") and '=' in arg:
                fp=open(arg.split('=')[1],'r')
                pdbA=fp.read()
                fp.close()
                pdb_list.append(pdbA)
            elif arg.startswith("-mol="):
                ter=arg[len("-mol="):]
            elif arg.startswith("-atom="):
                ter=arg[len("-atom="):]
            else:
                print("ERROR! No such option "+arg)

    if len(pdb_list)<2:
        print("ERROR! less than 2 chains")
        exit()

    sequence_list=[]
    for i in range(len(pdb_list)):
        sequenceA,pdbA=check_structure(pdb_list[i],atom=atom,mol=mol,
            chainID=chainID_list[i])
        sequence_list.append(sequenceA)
        pdb_list[i]=pdbA

    jobID=subprocess.Popen("date +%N",shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ).communicate()[0].strip()
    jobID="US"+jobID.decode()
   
    txt=''    
    for i in range(len(pdb_list)):
        txt+="%s%d\n"%(jobID,i+1)
        fp=open("%s/%s%d.pdb"%(outdir,jobID,i+1),'w')
        fp.write(pdb_list[i])
        fp.close()
    fp=open("%s/%s.list"%(outdir,jobID),'w')
    fp.write(txt)
    fp.close()
        
    print('<meta http-equiv="refresh" content="0.1;url=../output/%s.html">'%jobID)

    cmd=cmd_template.substitute(outdir=outdir,TMalign=TMalign,jobID=jobID,
        mol=mol,atom=atom)
    TMalign_output,stderr=subprocess.Popen(cmd,shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    TMalign_output=TMalign_output.decode()
    stderr=stderr.decode()

    make_display_structure(outdir, jobID)
    download_txt=''
    for i in range(len(pdb_list)):
        download_txt+=download_template.substitute(jobID=jobID, I=i+1)
    fp=open(outdir+'/'+jobID+".html",'w')
    fp.write(html_template.substitute(download_txt=download_txt,
        jobID=jobID, TMalign_output=TMalign_output.strip()))
    fp.write(stderr)
    fp.close()


    # failure to write bywhom will not affect user output display
    write_bywhom(jobID,'\t'.join([str(len(sequence)) for sequence in sequence_list]),
        '\t'.join(sequence_list))
