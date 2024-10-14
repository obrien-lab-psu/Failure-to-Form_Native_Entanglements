import glob
Finished_count = 0
pdbs = glob.glob('/storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Simulations_of_Native_Entanglement_Misfolding/Rebuild_AllAtom_structures/data/post_rebuilt/*')
for f in pdbs:
    tag = f.split('/')[-1].replace('_rebuilt.pdb', '')
    #print(f, tag)

    opt_path = f'../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Coarse_Graining_Model/{tag}/opt_nscal.log'
    #print(opt_path)

    # read in the opt_nscal.log and determine if there is a ## Final nscal: line indicating it finished
    opt_lines = [x for x in open(opt_path, 'r').readlines()]
    #print(opt_lines)

    Finished = False
  
    for line in opt_lines:
        if '## Final nscal:' in line:
            Finished = True
            Finished_count += 1
    print(f'{tag} Finished: {Finished}')
    
print(f'Finished_count: {Finished_count}')