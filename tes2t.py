def contact_check_multistep_v1(f_read, id_i, step1, step2, error_tolerence):

    df = pd.read_hdf(f_read, 'df')
    
    df_step = df.loc[df['step'].isin(list(range(step1, step2+1)))][[
        'step', 'id', 'type', 'x', 'y', 'z',
        'vx','vy','vz',
        'fx','fy','fz',
        'omegax','omegay','omegaz',
        'tqx','tqy','tqz',]]

    groups_byid = df_step.groupby(['id'])
    
    dfi = (groups_byid.get_group(id_i))
    
    typei = dfi[['type']].values[0:1]
    xi = dfi[['x','y','z']].values[:-1]
    vi = dfi[['vx','vy','vz']].values[:-1]
    fi = dfi[['fx','fy','fz']].values[:-1]
    omi = dfi[['omegax','omegay','omegaz']].values[:-1]
    tqi = dfi[['tqx','tqy','tqz']].values[:-1]

    atomi = single_atom(typei,xi,vi,fi,omi,tqi)

    print(xi)

    id_list = groups_byid.groups.keys()

    id_j_list = [i for i in id_list if i !=id_i]

    id_j_list = np.asarray(id_j_list)

    number_idj = len(id_j_list)
    number_wall = len(wall_list)
    # overlap array time in first dim. idj in second dim

    ifoverlap_ij_array = np.full((step2 - step1, number_idj), False)
    ifoverlap_iw_array = np.full((step2 - step1, number_wall), False)

    for n, id_j in enumerate(id_j_list):

        dfj = (groups_byid.get_group(id_j))
        
        def selectindexfori(dfj):
            indexjfori = (dfj[['step']].values[:-1]-step1).astype(int)
            indexjfori = indexjfori[:,0].T
            return indexjfori
        
        def selectstepfori(dfj):
            
            indexjfori = selectindexfori(dfj)
            typei = dfi[['type']].values[0:1]
            xi = dfi[['x','y','z']].values[indexjfori]
            vi = dfi[['vx','vy','vz']].values[indexjfori]
            fi = dfi[['fx','fy','fz']].values[indexjfori]
            omi = dfi[['omegax','omegay','omegaz']].values[indexjfori]
            tqi = dfi[['tqx','tqy','tqz']].values[indexjfori]
            atomi = single_atom(typei,xi,vi,fi,omi,tqi)
            
            return atomi

        typej = dfj[['type']].values[0:1]
        xj = dfj[['x','y','z']].values[:-1]
        vj = dfj[['vx','vy','vz']].values[:-1]
        fj = dfj[['fx','fy','fz']].values[:-1]
        omj = dfj[['omegax','omegay','omegaz']].values[:-1]
        tqj = dfj[['tqx','tqy','tqz']].values[:-1]

        atomj = single_atom(typej,xj,vj,fj,omj,tqj)

        jtoi = j_class(selectstepfori(dfj), atomj, 0, 0,0,1)
        
        

        ifoverlap_ij = jtoi.ifoverlap()
        overlapij_vector = jtoi.overlapij_vector()

        
        ifoverlap_ij_array[selectindexfori(dfj),n:n+1] = ifoverlap_ij

        
        print(ifoverlap_ij_array.dtype)
    
    id_j_wall_list = id_j_list

    for n, walllist in enumerate(wall_list):

        id_j_wall_list = np.append(id_j_wall_list, -n-1)
        if walllist[0] == 'p':
            wtoi = wall_p_class(
                atomi,walllist[1], walllist[2], 0,0,0,1,
                )
        elif walllist[0] == 'cy':
            wtoi = wall_cy_class(
                atomi,walllist[1], walllist[2], walllist[3], 0,0,0,1,
                )
        else:
            print('walltype not p not cy')
        print('sfwsfe')
        print(wtoi.atomi.x)

        ifoverlap_iw = wtoi.ifoverlap()
        print(ifoverlap_iw)
        print('kkk')
        print(wtoi.xij())
        print(wtoi.atomi.radius())
        print(wtoi.overlap_length())
        overlapiw_vector = wtoi.overlapij_vector()

        ifoverlap_iw_array[:,n:n+1] = ifoverlap_iw
        print(ifoverlap_iw_array.dtype)
    
    ifoverlap_ij_iw_array = np.concatenate((ifoverlap_ij_array, ifoverlap_iw_array), axis=1)
    print('ggg')
    print(ifoverlap_ij_iw_array)
    diff_next = (ifoverlap_ij_iw_array[:-1] != ifoverlap_ij_iw_array[1:])

    index_diff_next = np.nonzero(diff_next)
    print(index_diff_next)
    step_id_ifover_diffnext = np.empty((3, len(index_diff_next[0])), dtype=int)
    step_id_ifover_diffnext[2:3, :] = ifoverlap_ij_iw_array[index_diff_next[0], index_diff_next[1]]
    step_id_ifover_diffnext[0:1, :] = index_diff_next[0] + step1
    step_id_ifover_diffnext[1:2, :] = id_j_wall_list[index_diff_next[1]]
    step_id_ifover_diffnext = np.transpose(step_id_ifover_diffnext)

    return [step_id_ifover_diffnext, ifoverlap_ij_iw_array]