# Extra Credit Repository
##### First start by downloading the remote repository on github into the instance:
    cd labs
    git clone https://github.com/L5-BShi/Extra-Credit

## 1. Exploring BLAST and indentifying the homologs of our query gene

Download the query protein that will be used in BLASTp:

    ncbi-acc-download -F fasta -m protein XP_021368563.1

Make a directory titled AMPL in the Extra Credit Repository, this directory is where all the files for this research will be placed:

    mkdir /home/ec2-user/labs/Extra-Credit/AMPL

Next, we must put the query protein into BLASTp to obtain the homologs of this query protein:

    blastp -db ../allprotein.fas -query XP_021368563.1.fa -outfmt 0 -max_hsps 1 -out AMPL.blastp.typical.out
    
At the moment, the format is hard to read so we will perform another command to obtain a tabular output:

    blastp -db ../allprotein.fas -query XP_021368563.1.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle"  -max_hsps 1 -out AMPL.blastp.detail.out

There are far too many homologs resulting from the BLAST search, therefore to make the research more mangaeable while also decreasing the homologs present from chance, we will assign an e-value to filter the homologs:

    awk '{if ($6< 1e-35 )print $1 }' AMPL.blastp.detail.out > AMPL.blastp.detail.filtered.out

We are interested in the number of homologs in our filtered homologs list: That number is 41 found from this command below:

    wc -l AMPL.blastp.detail.filtered.out


We are also interested to see if each of the nine species has equal representation of data, therefore we used this command to determine the number of paralogs in each species:

    grep -o -E "^[A-Z][a-z]+\." AMPL.blastp.detail.filtered.out  | sort | uniq -c

**Paralogs in each Species:**
Aplanci 3
Avaga 9
Bbelicheri 5
Cintestinalis 2
Dmelanogaster 3
Egranulosus 2
Hsapiens 6
Lanatina 6
Myessoensis 5


## 2. To perform a multiple sequence alignment and to determine average percent identities of the homologs

Download alignbuddy from buddysite to later determine percent identity of the homologs including gaps

    pip install buddysuite
    sudo yum -y install gnu-free-mono-fonts.noarch

We will perform a global multiple sequence alignment in muscle using the filtered homologs as our input:

    muscle -in ~/labs/Extra-CreditA/AMPL/AMPL.homologs.fas -out ~/labs/Extra-Credit/AMPL/AMPL.homologs.al.fas

However, at the moment, the amino acid overlaps between homologs are hard to make out, therefore we use this command which color codes similar amino acids:

    alv -kli  ~/labs/lab4-$MYGIT/gqr/gqr.homologs.al.fas | less -RS

To further locate the more conserved amino acids, we use this function to only show 50% conserved amino acids in color:

    alv -kli --majority ~/labs/lab4-$MYGIT/gqr/gqr.homologs.al.fas | less -RS

In order to view this color-coded multiple sequence alignment in our remote repository, we must convert it to a pdf:

    muscle -in ~/labs/Extra-Credit/AMPL/AMPL.homologs.fas -html -out ~/labs/Extra-Credit/AMPL/AMPL.homologs.al.html
    sed -i 's/<PRE>/<pre style="font-family: 'FreeMono', monospaced;">/g' AMPL.homologs.al.html
    wkhtmltopdf ~/labs/Extra-Credit/AMPL/AMPL.homologs.al.html --print-media-type ~/labs/Extra-Credit/AMPL/AMPL.homologs.al.pdf

We will need to determine the width of the alignment with gaps included by using this command:

    alignbuddy  -al  ~/labs/Extra-Credit/AMPL/AMPL.homologs.al.fas

To determine the number of columns in the alignment that are invariant to discover a potential conensus sequence, we use this command:

    alignbuddy -dinv 'ambig' ~/labs/lab4-$MYGIT/AMPL/AMPL.homologs.al.fas | alignbuddy  -al

We used the t-coffee program to calculate the average percent identity of these homologs while excluding gaps to be 33.8%:

    t_coffee -other_pg seq_reformat -in ~/labs/lab4-$MYGIT/AMPL/AMPL.homologs.al.fas -output sim

We used the alignbuddy program to calculate the average percent identity of these homologs while including gaps to be 29.2%:

    alignbuddy -pi ~/labs/lab4-$MYGIT/AMPL/AMPL.homologs.al.fas | awk ' (NR>2)  { for (i=2;i<=NF  ;i++){ sum+=$i;num++} }
     END{ print(100*sum/num) } '

## 3. Determining the most optimal phylogeny using IQ-TREE and performing a mid-point rooting

In order to determine the most optimal phylogeny, we must use the program IQ-TREE. This program provide an unrooted tree:

    iqtree -s ~/labs/Extra-Credit/AMPL/AMPL.homologs.al.fas -bb 1000 -nt 2

We cannot view this unrooted tree unless we use this command. This command produces a Newick formatted tree with basal polytomy which is indicative of an unrooted tree:

    nw_display ~/labs/Extra-Credot/AMPL/AMPL.homologs.al.fas.treefile
    
We will use R to display the true unrooted tree without basal polytomy:

    Rscript --vanilla ~/labs/Extra-Credit/plotUnrooted.R  ~/labs/Extra-Credit/AMPL/AMPL.homologs.al.fas.treefile ~/labs/Extra-Credit/AMPL/AMPL.homologs.al.fas.treefile.pdf 0.4

**This is the unrooted tree:**
![](https://i.imgur.com/7Hb2y4S.png)
**It is very hard to comprehend as there are so many homologs yet so little space to put it.**

We are also interested in the mid-point rooted tree given by this command:

    gotree reroot midpoint -i ~/labs/Extra-Credit/AMPL/AMPL.homologs.al.fas.treefile -o ~/labs/Extra-Credit/AMPL/AMPL.homologs.al.mid.treefile

Similar to IQ-TREE, we cannot view this mid-point tree wihtout using a separate command, which will give the mid-point rooted tree in ASCII.

    nw_order -c n ~/labs/Extra-Credit/AMPL/AMPL.homologs.al.mid.treefile  | nw_display -

We want the mid-point rooted tree as a graphic in svg format instead of in ASCII, therefore we must use this command:

    nw_order -c n ~/labs/Extra-Credit/AMPL/AMPL.homologs.al.mid.treefile | nw_display -w 1000 -b 'opacity:0' -s  >  ~/labs/Extra-Credit/AMPL/AMPL.homologs.al.mid.treefile.svg -

**This is the Cladogram:**
![](https://i.imgur.com/pviQ4Ji.png)

A problem with the previous image is that the smaller branch lengths will be hard to see in the phylogram, therefore we can display the tree as a cladogram using this command:

    nw_order -c n ~/labs/Extra-Credit/AMPL/AMPL.homologs.al.mid.treefile | nw_topology - | nw_display -s  -w 1000 > ~/labs/Extra-Credit/AMPL/AMPL.homologs.al.midCl.treefile.svg -

## 4. Gene-species reconciliation trees using Notung and third kind to view gene duplication, speciation and loss events

In order to reconcile the gene and species tree, we must use the Notung program:

    java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -s ~/labs/Extra-Credit/species.tre -g ~/labs/Extra-Credit/AMPL/AMPL.homologs.al.mid.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/labs/Extra-Credit/AMPL/

However, an issue is that the internal nodes are not labelled, thereore we must use this command:

    grep NOTUNG-SPECIES-TREE ~/labs/Extra-Credit/AMPL/AMPL.homologs.al.mid.treefile.reconciled | sed -e "s/^\[&&NOTUNG-SPECIES-TREE//" -e "s/\]/;/" | nw_display -

**This is the gene-species reconciled tree using Notung:**
![](https://i.imgur.com/7J5ixi1.png)

We want to use third kind to place the gene-species reconcilation within the species tree so that the reconciliated tree is easier to comprehend, but third kind can only read xml, therefore we must convert the file to an RecPhyloXML object:

    python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/labs/Extra-Credit/AMPL/AMPL.homologs.al.mid.treefile.reconciled --include.species

Finally, we get to use third kind to create our gene-species reconciliated tree within the species tree:

    thirdkind -Iie -D 40 -f ~/labs/Extra-Credit/AMPL/AMPL.homologs.al.mid.treefile.reconciled.xml -o  ~/labs/Extra-Credit/AMPL/AMPL.homologs.al.mid.treefile.reconciled.svg

**This is the gene-species reconciliated tree within the species tree:**
![](https://i.imgur.com/4jP7mwy.png)

## 5. Predicting protein domains using RPS-BLAST to search the Pfam database

First, we must download the Pfam database in order to search it for domains that are present in our homologs:

    wget -O ~/data/Pfam_LE.tar.gz ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Pfam_LE.tar.gz && tar xfvz ~/data/Pfam_LE.tar.gz  -C ~/data

We will use RPS-BLAST in order to determine the protein domains that are in our protein homologs, then compare it to the Pfam database of conserved domains:

    rpsblast -query ~/labs/Extra-Credit/AMPL/AMPL.homologs.fas -db ~/data/Pfam -out ~/labs/Extra-Credit/AMPL/AMPL.rps-blast.out  -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .0000001

The next set of commands will plot the predicted Pfam domains on the phylogeny so now we should see the domains next to their respective homolog on the phylogeny. Furthermore, we will see the domains color coded and placed on a scale that shows the size of the domain and its placement within the protein to the total size of the protein:

    sudo /usr/local/bin/Rscript  --vanilla ~/labs/Extra-Credit/plotTreeAndDomains.r ~/labs/Extra-Credit/AMPL/AMPL.homologs.al.mid.treefile ~/labs/Extra-Credit/AMPL/AMPL.rps-blast.out ~/labs/Extra-Credit/AMPL/AMPL.homologs.fas   ~/labs/Extra-Credit/AMPL/AMPL.tree.rps.pdf
    
**This is the phylogentic tree with domains of each homolog to the right of them, furthermore, the different colors represent different domains:**

![](https://i.imgur.com/JJwyhVB.png)

