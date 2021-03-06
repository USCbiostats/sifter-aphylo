#!/usr/bin/python
import os
import pickle
import sys
import subprocess
import Queue
import threading
import getopt
import copy
import string
import glob
import csv
import re
import random

def usage():
    print "\n-----------------------------------------------------------------"
    print "Usage: "
    print "   sifter_run.py [options] <queries_prep_folder> <results_folder>"
    print "-----------------------------------------------------------------\n"
    print "Examples:"
    print "   sifter_run.py ../example/queries ../example/results\n"
    print "   sifter_run.py -n 10  ../example/queries ../example/results\n"
    print "   sifter_run.py -e PF12491,PF13820  ../example/queries ../example/results\n"
    print "   sifter_run.py -f PF12491,PF13820  ../example/queries ../example/results\n"
    print "   sifter_run.py -n 1 --if ../example/family_list.txt  ../example/queries ../example/results\n"
    print "This function runs SIFTER on the prepared files generated by 'sifter_prepare.py'."
    print "@author Sayed Mohammad Ebrahim Sahraeian (mohammad@compbio.berkeley.edu)"
    print "Please cite new paper:"
    print "-Sahraeian SME, Luo KR, Brenner SE (2015)"
    print "\nThe SIFTER algorithm presented in the following paper:"
    print "- Engelhardt BE, Jordan MI, Srouji JR, Brenner SE. 2011. Genome-scale phylogenetic function annotation of large and diverse protein families. Genome Research 21:1969-1980. \n"
    print "inputs:"
    print "        <queries_prep_folder>    Path to the queries preperation."
    print "                                 folder. The folder should contain the" 
    print "                                 necessary query files (generated"
    print "                                 by 'sifter_prepare.py' script). Use" 
    print "                                 the same directory used as output"
    print "                                 in 'sifter_prepare.py' script" 
    print "        <results_folder>         Path to the output folder where"
    print "                                 results will be written to." 
    print "options:"
    print "           -n          INT       Number of threads (Default=4)"
    print "           -e          STRING    List of families for which you want"
    print "                                 to exclude running SIFTER on."
    print "                                 (in comma seperated format)"
    print "           --ie        STRING    Path to the input file where the"
    print "                                 list of families for which you"
    print "                                 want to exclude running SIFTER"
    print "                                 on."
    print "           -f          STRING    List of families for which you want"
    print "                                 to run SIFTER on."
    print "                                 (in comma seperated format)"
    print "                                 If not provided, SIFTER will run on"
    print "                                 all families in queries_prep_folder"
    print "           --if        STRING    Path to the input file where the"
    print "                                 list of families for which you"
    print "                                 want to run SIFTER on."
    print "                                 If not provided, SIFTER will run on"
    print "                                 all families in queries_prep_folder"
    print "           -t          INT       Number of functions to truncate" 
    print "                                 to in approximation [Default:"
    print "                                 Automatically computed in sifter_prepare.py]"
    print "                                 Smaller value leads to faster running time."
    print "           -h                    Help. Print Usage."



def write_goa_anns_to_pli_constrained(evidence_file,goa_anns, fam_id, seq_lookup,evidence_constraints):
    '''
    This converts the database rows to B. Engelhardt's arbitrary evidence XML foramt.
    Input looks like:
            {'A2VE79': [{'acc': 'GO:0000287',
             'code': 'ISS',
             'full_name': 'Diphosphoinositol polyphosphate phosphohydrolase 1',
             'genus': 'Bos',
             'is_not': 0L,
             'name': 'magnesium ion binding',
             'species': 'taurus',
             'symbol': 'NUDT3',
             'xref_dbname': 'UniProtKB',
             'xref_key': 'A2VE79'},
            {'acc': 'GO:0008486',
             'code': 'ISS',
             'full_name': 'Diphosphoinositol polyphosphate phosphohydrolase 1',
             'genus': 'Bos',
             'is_not': 0L,
             'name': 'diphosphoinositol-polyphosphate diphosphatase activity',
             'species': 'taurus',
             'symbol': 'NUDT3',
             'xref_dbname': 'UniProtKB',
             'xref_key': 'A2VE79'},
            ...
    '''
    f = open(evidence_file, 'w')
    f.write("<?xml version=\"1.0\"?>\n<Family>\n")
    f.write("  <FamilyID>%s</FamilyID>\n"%fam_id)
    
    for p_id, anns in goa_anns.iteritems():
        filtered_anns=[w for w in anns if w['code'] in evidence_constraints]
        if not filtered_anns:
            continue
        f.write("  <Protein>\n")
        f.write("    <ProteinName>%s</ProteinName>\n"%seq_lookup[p_id])
        f.write("    <ProteinNumber>%s</ProteinNumber>\n"%p_id)
        go_str = ''
        moc_str = ''
        for i,a in enumerate(filtered_anns):
            go_str += a['acc'][3:]
            moc_str += a['code']
            if i < len(filtered_anns)-1:
                go_str += ', '
                moc_str += ', '
        f.write("    <GONumber>%s</GONumber>\n"%('['+go_str+']'))
        f.write("    <MOC>%s</MOC>\n"%('['+moc_str+']'))
        f.write("  </Protein>\n")
    f.write("</Family>\n")
    f.close()
    
class ProcessingThread(threading.Thread):
    """Thread for running sequence alignments on a given input homolog cluster."""
    def __init__(self, thread_queue):
        threading.Thread.__init__(self)
        self.thread_queue = thread_queue
    def thread_operation(self, thread_data):
        query_data, output_prefix,pfam_id = thread_data
        self.output_prefix = output_prefix
        try:
            print "--------------------------------------------------"
            print "Running SIFTER for ",query_data['pfam_id'],':',pfam_id
            
            # Input evidence
            rand_id_1=random.randint(1000000,9999999)                    
            rand_id_2=random.randint(1000000,9999999)                    
            evidence_file = query_data['annotation_loc']
            evidence_file_pickle = query_data['annotation_loc_pickle']
            evidence_format = query_data['annotation_format']
            evidence_constraints = query_data['evidence_constraints']
            tree_file = query_data['tree_loc']
            if os.path.exists(tree_file+'.gz'):
                if os.path.exists('%s.%d'%(tree_file,rand_id_1)):
                    subprocess.check_call("rm %s.%d"%(tree_file,rand_id_1),shell=True)                    
                subprocess.check_call("gunzip -c %s.gz > %s.%d"%(tree_file,tree_file,rand_id_1),shell=True)
            if os.path.exists(evidence_file_pickle+'.gz'):
                if os.path.exists('%s.%d'%(evidence_file,rand_id_2)):                
                    subprocess.check_call("rm %s.%d"%(evidence_file,rand_id_2),shell=True)                    
                if os.path.exists('%s.%d'%(evidence_file_pickle,rand_id_2)):                
                    subprocess.check_call("rm %s.%d"%(evidence_file_pickle,rand_id_2),shell=True)                    
                subprocess.check_call("gunzip -c %s.gz > %s.%d"%(evidence_file_pickle,evidence_file_pickle,rand_id_2),shell=True)
                [ee,anns, p, seq_lookup]=pickle.load(open("%s.%d"%(evidence_file_pickle,rand_id_2), 'r'))
                write_goa_anns_to_pli_constrained("%s.%d"%(evidence_file,rand_id_2),anns, p, seq_lookup,evidence_constraints)
            pfam_id = query_data['pfam_id']
            n_terms=query_data['n_terms']
            if query_data['enforced_trunc']>0:
                max_simul_fcns=query_data['enforced_trunc']
            else:
                max_simul_fcns=query_data['e_time'][0]
                
            print "n_terms :",n_terms
            print "truncation level :",max_simul_fcns
            print "Estimated running time for family %s = %s (95%% confidence upper bound = %s)"%(pfam_id,query_data['e_time'][4],query_data['e_time'][5])            
            evidence_str = ''
            for ev_type in evidence_constraints:
                evidence_str = evidence_str + "--with-" + string.lower(ev_type) + " "
            go_ontology_sqlite = main_dir+'/data/goterms.sqlite'



            java_sifter_cmd = \
                "java -jar -Xmx4g " + sifter_java + " " +\
                "--generate " + " " +\
                evidence_str + " " +\
                "--protein " + evidence_file+'.%d'%rand_id_2+ " " +\
                "--reconciled " + tree_file+'.%d'%rand_id_1 + " " +\
                "--ontology " + go_ontology_sqlite + " " +\
                "--output " + self.output_prefix+'_result.txt ' +\
                pfam_id
            # Copy xml into the "input directory" of the new evaluation directory
            print "[aphylo] Executing the following Command:\n"
            print java_sifter_cmd;
            exec_command = java_sifter_cmd
            #print "Running: "+exec_command
            retcode = subprocess.call(exec_command, shell=True)
            if (retcode == 0):
                print "Success: Ran generation step.\n"
            else:
                print "Error: Sifter failed to generate.\n"
                raise Exception
                    
          
            java_sifter_cmd = \
                "java -jar -Xmx4g " + sifter_java + " " +\
                evidence_str + " " +\
                "--protein " + evidence_file+'.%d'%rand_id_2 + " " +\
                "--reconciled " + tree_file+'.%d'%rand_id_1 + " " +\
                "--ontology " + go_ontology_sqlite + " " +\
                "--output " + self.output_prefix +'_result.txt ' +\
                "--familyfile " + self.output_prefix[:-len(pfam_id)] + "/infer-" + pfam_id + ".fx " +\
                "--scale " + self.output_prefix[:-len(pfam_id)] + "/scale-" + pfam_id + ".fx " +\
                "--alpha " + self.output_prefix[:-len(pfam_id)] + "/alpha-" + pfam_id + ".fx " +\
                ("--truncation %d "%max_simul_fcns) +\
                "--xvalidation --folds 0 " +\
                pfam_id
           
            print "[aphylo] Executing the following command:\n" 
            print java_sifter_cmd
            exec_command = java_sifter_cmd
            #print "Running: "+exec_command
            retcode = subprocess.call(exec_command, shell=True)
            if (retcode == 0):
                print "Success: Ran inference step.\n"
            else:
                print "Error: Sifter failed to run inference.\n"            


            res_file=self.output_prefix+'_result.txt'
            if os.path.exists(res_file):    
                results={}
                data = list(csv.reader(open(res_file, 'rb'), delimiter='\t')) 
                if data:
                    candids = []
                    store_data = {}
                    candids=data[0][1:]
                    for line in data[1:]:
                        gene=line[0].split('/')[0]
                        if gene not in results:
                            results[gene]=[]
                        if len(line[0].split('/'))>1:
                            results[gene].append([line[0].split('/')[1],[float(w) for w in line[1:]]])
                        else:
                            results[gene].append(['',[float(w) for w in line[1:]]])  
                    res_pickle_file=self.output_prefix+'_result.pickle'
                    pickle.dump({'results':results,'candids':candids},open(res_pickle_file,'w')) 
            
            
            if os.path.exists('%s.%d'%(evidence_file,rand_id_2)):                
                subprocess.check_call("rm %s.%d"%(evidence_file,rand_id_2),shell=True)                    
            if os.path.exists('%s.%d'%(evidence_file_pickle,rand_id_2)):                
                subprocess.check_call("rm %s.%d"%(evidence_file_pickle,rand_id_2),shell=True)                    
            if os.path.exists('%s.%d'%(tree_file,rand_id_1)):                
                subprocess.check_call("rm %s.%d"%(tree_file,rand_id_1),shell=True)                    
            
            
            print "Query complete for ", pfam_id
        
        except Exception as e:
            print >> sys.stderr, "Error evaluating %s"%query_data['pfam_id']
            print >> sys.stderr, "Error: ", e
            exit(1)
        
    def flag_start(self):
        f = open(self.output_prefix + ".sifterj.processing", 'w')
        f.close()
    
    def unflag_start(self):
        os.remove(self.output_prefix + ".sifterj.processing")
    
    def flag_finish(self):
        #self.unflag_start()
        f = open(self.output_prefix + ".sifterj.processed", 'w')
        f.close()
    
    def run(self):
        while True:
            # Spawn a thread with data from the queue
            thread_data = self.thread_queue.get()
            self.output_prefix = thread_data[1]
            
            # Run thread's function on the data
            try:
                self.flag_start()
                self.thread_operation(thread_data)
                self.flag_finish()
            except:
                self.unflag_start()
                print "Unexpected thread error:", sys.exc_info()[0]
                print "Thread data:", thread_data
            # Send signal that this task finished
            self.thread_queue.task_done()


    
if __name__=="__main__":
    
    # Initialization
    main_dir=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    sifter_java=os.path.dirname(main_dir)+'/core/sifter2.1.1.jar'
    if not os.path.exists(sifter_java):
        print "\nERROR: No SIFTER jar file exists at %s\n"%sifter_java
        sys.exit()
        
    num_threads=4
    exclude_fams=[]
    # Check for options
    
    exclude_fams=[]
    only_fams=[]
    opts, args = getopt.getopt(sys.argv[1:], "hn:e:f:t:",['ie=','if=']) 
    truncation_level=0
    
    if len(args) != 2:
        usage()
        sys.exit()
    if len(opts)>0:
        for o, a in opts:
            if o == "-n":
                num_threads=int(a)
            elif o == "-e":
                splited =a.strip().split(',')
                exclude_fams.extend(list(set([w for w in splited if w])))
            elif o == "--ie":
                x_input_file = a        
                if not os.path.exists(x_input_file):
                    print "\nERROR: No file exists at %s\n"%x_input_file
                    sys.exit()
                f = open(x_input_file, 'r')
                a=f.read()
                splited =re.split(' |,|;|\n',a.strip())
                exclude_fams.extend(list(set([w for w in splited if w])))                
            elif o == "-f":
                splited =a.strip().split(',')
                only_fams.extend(list(set([w for w in splited if w])))
            elif o == "--if":
                i_input_file = a        
                if not os.path.exists(i_input_file):
                    print "\nERROR: No file exists at %s\n"%i_input_file
                    sys.exit()
                f = open(i_input_file, 'r')
                a=f.read()
                splited =re.split(' |,|;|\n',a.strip())
                only_fams.extend(list(set([w for w in splited if w])))                
            elif o == "-t":
                truncation_level=int(a)
            else:
                usage()
                sys.exit()
   
    queries_folder=args[0]   
    if not os.path.exists(queries_folder):
        print "\nERROR: Queries preperation folder ( %s ) does not exist\n"%queries_folder
        print "Please use the same directory used as output in 'sifter_prepare.py'"
        sys.exit()
    results_folder=args[1]   
    if not os.path.exists(results_folder):
        os.mkdir(results_folder)
        
    print "\n-----------Reading the queries information-------------\n"
    query_files=glob.glob(queries_folder+"/*_query.pickle")
    queries_to_process = []   
    for qfile in query_files:
        pfam_id=qfile.split('/')[-1].split('_')[0]
        if only_fams:
            if pfam_id not in only_fams:
                continue                
        if pfam_id in exclude_fams:
            print "Excludeing family %s"%pfam_id
            continue
        output_prefix = results_folder + '/' + pfam_id    
        if  not(os.path.isfile(output_prefix + "_result.pickle")):
            query_data = pickle.load(open(qfile, "rb" ))
            query_data['output_to'] = output_prefix
            query_data['enforced_trunc']=truncation_level            
            queries_to_process.append((copy.deepcopy(query_data), output_prefix,pfam_id))
        else:
            print "Predictions have already been made for %s"%pfam_id
            
    print "\n-------------------Running SIFTER-----------------------\n"
    thread_queue = Queue.Queue()
    for i in range(num_threads):
        t = ProcessingThread(thread_queue)
        t.setDaemon(True)
        t.start()
     
    for thread_data in queries_to_process:
        thread_queue.put(item=thread_data, block=False)
    
    # Wait on the queue until everything has been processed         
    thread_queue.join()


    errors=0
    qs=0
    all_qs=[]
    for w in queries_to_process:        
        pfam_id=w[2]
        all_qs.append(pfam_id)
        if pfam_id in exclude_fams:
            continue        
        if only_fams:
            if pfam_id not in only_fams:
                continue                
        qs+=1
        res_pickle_file=w[1]+'_result.pickle'
        if not (os.path.isfile(res_pickle_file)):
            errors+=1

    if only_fams:
        all_qs=only_fams
    
    all_qs=set(all_qs)-set(exclude_fams)
    n_qs=len(all_qs)-qs
    succus=qs-errors

    print "\nSIFTER results are ready for %d out of %d families. (%s missed due to errors, %s missed due to lack of query files)"%(succus,len(all_qs),n_qs,errors)

    if succus>0:  
        print "-------------------Runnig SIFTER is Done----------------------"
        print "You may extract results for specific queries using  'sifter_extract.py'.\n"

