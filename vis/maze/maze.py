from __future__ import print_function
import re, os, glob, sys
from tempfile import NamedTemporaryFile
import itertools as itt
from subprocess import check_output, call, Popen, PIPE


################################################################################
#
# CONFIG
#
# Specify the paths to external tools, if not in PATH
#
tool = {}
tool['mummer']     = 'mummer'
tool['lastdb']     = 'lastdb' # last >= 584 is needed. Todo(meiers): write requirements
tool['lastal']     = 'lastal'
tool['last-split'] = 'last-split'
#
# END CONFIG
#
################################################################################


m_pat = re.compile(r'\s*(\d+)\s+(\d+)\s+(\d+)')

def mummer_matches(fn1, fn2, l, m='mumreference', verbose=False):
    cmd = '{} -{} -l {} -b -c {} {} 2>/dev/null'.format(tool['mummer'], m, l, fn1, fn2)
    if verbose:
        print('#', cmd)
    mummer_out = check_output(cmd, shell=True)
    fwd = []
    rev = []
    for l in mummer_out.splitlines():
        if l.startswith('>'):
            is_rev = True if 'Reverse' in l else False
        else:
            s1, s2, l = [int(x) for x in m_pat.match(l).group(1,2,3)]
            if is_rev:
                rev.append([s1, s1+l-1, s2, s2-l+1])
            else:
                fwd.append([s1, s1+l-1, s2, s2+l-1])
    return {'fwd': fwd, 'rev': rev}


################################################################################
# maze detail 

def LASTsplit_matches(ref, seq):
    '''Given two DNA sequences, compute last-split alignment
    and return matches'''
    matches = []
    with NamedTemporaryFile(delete=False, suffix='.fa') as f_ref, \
         NamedTemporaryFile(delete=False, suffix='.fa') as f_read:
        fn_ref = f_ref.name
        fn_read = f_read.name
        print('>tmp_ref', file=f_ref)
        print(ref, file=f_ref)
        print('>tmp_read', file=f_read)
        print(seq, file=f_read)
    try: 
        cmd = '{} {} {} 2>/dev/null'.format(tool['lastdb'], fn_ref, fn_ref)
        call(cmd, shell=True)
        p1 = Popen([tool['lastal'], '-a', '3', fn_ref, fn_read], stdout=PIPE)
        p2 = Popen(tool['last-split'], stdin=p1.stdout, stdout=PIPE)
        last_split = p2.communicate()[0]
    except Exception as e:
        print ("FATAL: LAST crashed.", e, sep='\n', file=sys.stderr)
        last_split=''
    finally: 
        #for x in glob.glob(str(fn_ref + '*')):
        #    os.remove(x)
        #os.remove(fn_read)
        pass
    matches = []
    it = iter(last_split.splitlines())
    while True:
        try:
            ln = it.next()
        except StopIteration, e:
            break
        if ln.startswith("s"):
            ln2 = it.next()
            match=dict()
            # reference entry
            s = ln.split()
            match['d1'] = int(s[2])
            match['d2'] = int(s[2])+int(s[3])
            match['ali1'] = s[6]
            # query entry
            s = ln2.split()
            match['q1'] = int(s[2])
            match['q2'] = int(s[2])+int(s[3])
            match['ali2'] = s[6]
            match['strand'] = s[4]
            match['qlen'] = int(s[5])
            # mismatches
            match['sim'] = len(match['ali1']) - sum([a==b for (a,b) in 
                            zip(match['ali1'].upper(), match['ali2'].upper())])
            match['record'] = _format_last_alignment(match)
            matches.append(match)

    return matches

def _format_last_alignment(m, w=60):
    a1 = m['ali1']
    a2 = m['ali2']
    am = ''.join([ (':' if c.upper()==C.upper() else ' ') for c,C in zip(a1,a2) ])
    l1 = []
    l2 = []
    lm = []
    while a1:
        l1.append(a1[:w])
        l2.append(a2[:w])
        lm.append(am[:w])
        a1 = a1[w:]
        a2 = a2[w:]
        am = am[w:]
    final = ''        
    for i in range(len(l1)):
        final += 'Ref   +'.format(m['strand']).ljust(10) + l1[i] + '\n'
        final += ''.ljust(10) + lm[i] + '\n'
        final += 'Query {}'.format(m['strand']).ljust(10) + l2[i] + '\n'
        final += '\n'
    return (final)

def _transform_coords(match):
    '''Transform coordinates from LAST type (0-based, minus strand are 
        coordinates in the reverse complement sequence) to traditional
        coordinates (1-based, end is last based covered)'''
    m = match.copy()
    m['d1'] = match['d1'] + 1
    if m['strand'] == '+':
        m['q1'] = match['q1'] + 1
    else:
        m['q1'] = match['qlen'] - match['q2'] + 1
        m['q2'] = match['qlen'] - match['q1']
    return m

def LASTsplit_breakpoints(ref, seq, matches, W=40, V=20): # matches must be sorted by query
    breakpoints = []
    for m1,m2 in zip(matches, matches[1:]):
        just = 16
        print_qu_coords =   (m1['q2'] if m1['strand']=='+' else len(seq) - m1['q1'],\
                            (m2['q1'] if m2['strand']=='+' else len(seq) - m2['q2'])+1)

        S    = ref      if m1['strand'] == '+' else _rc(ref)
        T    = ref      if m2['strand'] == '+' else _rc(ref)
        a    = m1['d2'] if m1['strand'] == '+' else len(ref) - m1['d1']
        b = m1['q2'] if m1['strand']=='+' else len(seq) - m1['q1']
        c = m2['q1'] if m2['strand']=='+' else len(seq) - m2['q2']
        d    = m2['d1'] if m2['strand'] == '+' else len(ref) - m2['d2']

        # upper panel
        a_pr = m1['d2'] if m1['strand'] == '+' else m1['d1']
        upper ="Ref{} {}-{}:".format(m1['strand'], 
                                     a_pr-W+1, 
                                     a_pr+V).ljust(just) +\
                S[a-W : a].upper().rjust(W) +\
                S[a   : a+V].lower()

        # middle panel
        spacer = seq[b   : c  ].lower() if c-b < 20 else (seq[b:b+5] + '...' + seq[c-5:c]).lower()
        middle = "Read {}-{}:".format(b-W+1, c+W).ljust(just) +\
                seq[b-W : b  ].upper().rjust(W) +\
                spacer.lower() +\
                seq[c   : c+W].upper()
        
        # lower panel:
        d_pr = m2['d1'] if m2['strand'] == '+' else m2['d2']
        lower = "Ref{} {}-{}:".format(m2['strand'],
                                      d_pr-V+1,
                                      d_pr+W).ljust(just) +\
                (T[d-V : d  ].lower() +\
                 T[d   : d+W].upper()).rjust(2*W+len(spacer))

        # HARDCODED MIHO LIMIT: 50
        mm_left = _miho(T[max(d-50,0):d][::-1], seq[max(c-50,0):c][::-1])
        mm_right = _miho(S[a:a+25], seq[b:b+25])

        # HTML outuput
        is_miho = print_qu_coords[1] - print_qu_coords[0] <= 1
        html  = upper[:just] + '<span><b>' + upper[just:just+W] + '</b>' + upper[just+W:just+W+mm_right] + '</span>' + upper[just+W+mm_right:] + '\n'
        if is_miho:
            html += middle[:just] + '<span><b>' + middle[just:] + '</b></span>\n'
        else:
            html += middle[:just] + '<span><b>' + middle[just:just+W] + '</b></span>' + middle[just+W:just+W+len(spacer)] + '<span><b>' + middle[just+W+len(spacer):] + '</b></span>\n'
        html += lower[:max(just,len(lower)-W-mm_left)] + '<span>' + lower[max(just,len(lower)-W-mm_left):max(just,len(lower)-W)] + '<b>' + lower[max(just,len(lower)-W):] + '</span>' + '\n'
        
        breakpoints.append(dict(html=html,
                                qu_coords=(print_qu_coords[0], print_qu_coords[1]),
                                match_index=str(matches.index(m1))+','+str(matches.index(m2)), 
                                event='Microhomology' if is_miho else 'Insertion',
                                length=mm_left + mm_right if is_miho else print_qu_coords[1] - print_qu_coords[0] - 1))
    return breakpoints

def _rc(seq):
    rc = dict(a='t', c='g', g='c', t='a',
              A='T', C='G', G='C', T='A')
    return "".join([rc[b] if b in rc else 'N' for b in seq[::-1]])

def _miho(seq1, seq2):
    count=0
    for x in [a==b for (a,b) in zip(seq1, seq2)]:
        if not x:
            break
        count += 1
    return count


def LASTsplit_ref_breakpoints(ref, seq, matches, W=40, V=25):
    X=100
    breakpoints = []
    for m1, m2 in itt.combinations(matches, 2):
        # m2 is always >= m1 on the query
        if abs(m1['d2'] - m2['d1']) < X:
            # exclude matches with sequencing end
            if matches.index(m1) == 0 and m1['strand']=='-':
                continue
            left, right = m1, m2
        elif abs(m2['d2'] - m1['d1']) < X:
            # exclude matches with sequencing end
            if matches.index(m2) == len(matches)-1 and m2['strand']=='+':
                continue
            left , right = m2, m1
        else:
            continue
        just = 16
        S = seq if left['strand']  == '+' else _rc(seq)
        T = seq if right['strand'] == '+' else _rc(seq)
        a = left['q2']
        b = left['d2']
        c = right['d1']
        d = right['q1']

        # upper panel
        pr_a = left['q2'] if left['strand'] == '+' else len(seq)-left['q2']
        upper ="Qu{} {}-{}:".format(left['strand'], pr_a-W+1, pr_a+V).ljust(just) +\
                S[a-W : a].upper().rjust(W) +\
                S[a   : a+V].lower()

        # middle panel
        spacer = ref[b   : c  ] #if c-b < 20 else seq[b:b+5] + '...' + seq[c-5:c]
        middle = "Ref {}-{}:".format(b-W+1, c+W).ljust(just) +\
            ref[b-W : b  ].upper().rjust(W) +\
            spacer.lower() +\
            ref[max(b,c)   : max(b,c)+W].upper().ljust(W) # at least length 2*W

        # lower panel
        pr_d = right['q1'] if right['strand'] == '+' else (len(seq)-right['q1'])
        lower = "Qu{} {}-{}:".format(right['strand'], pr_d-V+1, pr_d+W).ljust(just) +\
                (T[d-min(V, W-b+c) : d ].lower() +\
                T[d   : d+W+max(b-c,0)].upper()).rjust(len(middle)-just)

        # HARDCODED MIHO LIMIT: 50 (!)
        mm_right = _miho(S[a:a+50], ref[b:b+50])
        mm_left = _miho(T[d-50:d][::-1], ref[c-50:c][::-1])

        # HTML outuput
        is_overlap = b>c
        html  = 'Microhomology:' + str((mm_right,mm_left)) + '\n'
        html += upper[:just] + '<span><b>' + upper[just:just+W] + '</b>' + upper[just+W:just+W+mm_right] + '</span>' + upper[just+W+mm_right:] + '\n'

        # middle
        html         += middle[:just]
        html         += '<span><b>'
        if is_overlap:
            html     += middle[just:]
        else:
            html     += middle[just : just+W]
            html     += '</b>' 
            # no consecutive coloring
            if mm_left + mm_right < len(spacer):
                html += middle[just+W          : just+W+mm_right]
                html += '</span>'
                html += middle[just+W+mm_right : just+W+len(spacer)-mm_left]
                html += '<span>'
                html += middle[just+W+len(spacer)-mm_left : just+W+len(spacer)]
            # consecutive coloring
            else: 
                html += middle[just+W : just+len(spacer)]
            html     += '<b>'
            html     += middle[just+W+len(spacer) : ]
        html     += '</b></span>' + '\n'
        html         += lower[:max(just,len(lower)-W-mm_left-max(0,b-c))] + '<span>' + lower[max(just,len(lower)-W-mm_left-max(0,b-c)):max(just,len(lower)-W-mm_left-max(0,b-c))] + '<b>' + lower[max(just,len(lower)-W-mm_left-max(0,b-c)):] + '</span>' + '\n'

        breakpoints.append(dict(html=html,
                                match_index=str(matches.index(m1))+','+str(matches.index(m2)), 
                                ref_coords=(b+1,c+1),
                                event='Overlap' if is_overlap else 'Distance',
                                miho=(mm_left, mm_right),
                                length=abs(c-b)))
    return breakpoints

