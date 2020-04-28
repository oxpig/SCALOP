# SCALOP

DESCRIPTION
 
    SCALOP
    Sequence-based antibody Canonical LOoP structure annotation
 
    o To setup paths etc on your machine run:
        python setup.py install
    
    o Input format:
        o Single complete antibody sequence 
        o Complete heavy and light chain sequences separated by '/'
        o Fasta file 

    o Python API
        o The SCALOP python package can be imported using:
        >>> from scalop.predict import assign
        o Results are shown as follows:
        {'input':                                              # Input sequence name 
            {'H1':                                             # CDR type
                ['H1', 'TASGFNIKDYYIH', 'H1-13-A', '4hie_B'],  # CDR type, CDR on the input sequence, assigned cluster, median structure of the cluster
             'H2': 
                ['H2', 'WIDPEIGDTE', 'H2-10-A', '1aj7_H']}
        }
        
    o Command line 
        o To run SCALOP on a command line:
        $ SCALOP -i <input_sequence> --scheme <numbering_scheme> --definition <cdr_definition>
        
    o Pre-requisites:
        o HMMER 3.1b2 (February 2015); http://hmmer.org/
        o numpy v. 1.10+
        o pandas 
		
AUTHORS
 
	2018
	Wing Ki Wong
	Dr Jinwoo Leem
	Prof Charlotte M. Deane - Oxford Protein Informatics group.
	
	Contact: opig@stats.ox.ac.uk


