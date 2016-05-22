import itertools
import os
import numpy as np

AMINO_ACID_LIST = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
                   'T', 'V', 'W', 'Y', 'X']
AMINO_ACID_CODES = {code: i + 1 for (i, code) in enumerate(AMINO_ACID_LIST)}
TOTAL_TRIPLETS = [''.join(tup) for tup in itertools.combinations(AMINO_ACID_LIST[:-1], 3)]

def run_haystack():
    print 'Loading Patient Data'
    patients_data = load_patients_data()
    dom_frac = .6

    dom_thresh = int(len(patients_data['MS']) * dom_frac)

    filtered_triplets = []
    atomic_motifs = []

    print 'Generating Atomic Motifs'
    for triplet in TOTAL_TRIPLETS:
        motif = maximize_atomic_motif(patients_data, triplet)
        score_dict = score_motif_match(patients_data, motif)
        if score_dict['ms_count'] >= dom_thresh:
            filtered_triplets.append(triplet)
            atomic_motifs.append(motif)

    print "Generated %s total Atomic Motifs" % len(atomic_motifs)

    total_doms = []

    for i, motif in enumerate(atomic_motifs):
        print 'Running local search from Atomic Motif %d'%i
        print motif
        discovered_doms = local_search(motif, patients_data, filtered_triplets,
                                       dom_thresh=dom_thresh,
                                       atomic_motifs=atomic_motifs)
        if not discovered_doms:
            print 'No DOMs discovered'
        else:
            total_doms.extend(discovered_doms)
            print '%d new Disease-Only-Motifs of %d total DOMs discovered'%(len(discovered_doms),
                                                                            len(total_doms))
    return total_doms

def local_search(starting_motif, patients_data, filtered_triplets, gap_count=1,
                 dom_thresh=27, atomic_motifs=None):

    score = score_motif_match(patients_data, starting_motif)
    if score['ms_count'] >= dom_thresh:
        if score['hc_count'] == 0:
            print 'Discovered DOM'
            print "Matches %d MS patients"%(score['ms_count'])
            return [(starting_motif, score['ms_count'])]

    queue = [starting_motif]
    while queue:
        motif = queue.pop(0)
        # Executing EMS Extension
        doms, new_extensions = get_ems_extensions(patients_data, motif, filtered_triplets,
                                                  dom_thresh=dom_thresh)

        if doms:
            return doms
        if new_extensions:
            queue.extend(new_extensions)
        if len(motif) <= gap_count:
            # Executing Cardinality Extension
            doms, new_extensions = get_cardinality_extensions(patients_data, motif,
                                                              filtered_triplets,
                                                              atomic_motifs=atomic_motifs,
                                                              dom_thresh=dom_thresh)

            if doms:
                return doms
            if new_extensions:
                queue.extend(new_extensions)
    return []


def maximize_atomic_motif(patients_data, triplet, start_range=-100, end_range=10):

    motifs = [(triplet, i) for i in xrange(start_range, end_range)]
    scores = [score_motif_match(patients_data, [e]) for e in motifs]
    best_score = max(scores, key=lambda x: x['ms_count'])
    return [(triplet, range(start_range, end_range)[scores.index(best_score)])]

def get_cardinality_extensions(patients_data, motif, filtered_triplets, atomic_motifs=None,
                               start_range=-100, end_range=10, dom_thresh=27, score_thresh=15):
    bad_start_indices = []
    for (explicit_motif_subseq, rel_position) in motif:
        bad_start_indices.extend(range(rel_position, rel_position + len(explicit_motif_subseq)))

    discovered_doms = []
    extended_motifs = []
    if atomic_motifs is None:
        atomic_motifs = []
        for triplet in filtered_triplets:
            atomic_motifs.extend([[(triplet, i)] for i in xrange(start_range, end_range)
                                  if i not in bad_start_indices])

    for atomic_motif in atomic_motifs:
        extended_motif = motif + atomic_motif
        score = score_motif_match(patients_data, extended_motif)
        if score['score'] >= score_thresh and score['ms_count'] >= dom_thresh:
            if score['hc_count'] == 0:
                print 'Discovered DOM'
                print extended_motif
                print "Matches %d MS patients"%(score['ms_count'])
                discovered_doms.append((extended_motif, score['ms_count']))
            else:
                extended_motifs.append(extended_motif)

    return (discovered_doms, extended_motifs)


def get_ems_extensions(patients_data, motif, filtered_triplets, start_range=-100, end_range=10,
                       dom_thresh=27, score_thresh=15):
    bad_start_indices = []
    for (explicit_motif_subseq, rel_position) in motif:
        bad_start_indices.extend(range(rel_position, rel_position + len(explicit_motif_subseq)))
    new_motifs = []

    for i, (explicit_motif_subseq, rel_position) in enumerate(motif):
        back_position = rel_position - 1
        if back_position < start_range or back_position in bad_start_indices:
            continue
        for aa in AMINO_ACID_LIST:
            new_ems = aa + explicit_motif_subseq
            if new_ems[:3] not in filtered_triplets:
                continue
            extended_motif = motif[:]
            extended_motif[i] = (new_ems, back_position)
            new_motifs.append(extended_motif)

        forward_position = rel_position + 1
        if forward_position >= end_range or forward_position in bad_start_indices:
            continue

        for aa in AMINO_ACID_LIST:
            new_ems = explicit_motif_subseq + aa
            if new_ems[-3:] not in filtered_triplets:
                continue
            extended_motif = motif[:]
            # Do not update `rel_position` extending for forward extension
            extended_motif[i] = (new_ems, rel_position)
            new_motifs.append(extended_motif)

    discovered_doms = []
    extended_motifs = []
    for extended_motif in new_motifs:
        score = score_motif_match(patients_data, extended_motif)
        if score['score'] >= score_thresh and score['ms_count'] >= dom_thresh:
            if score['hc_count'] == 0:
                print 'Discovered DOM'
                print extended_motif
                print "Matches %d MS patients"%(score['ms_count'])
                discovered_doms.append((extended_motif, score['ms_count']))
            else:
                extended_motifs.append(extended_motif)

    return (discovered_doms, extended_motifs)


def score_motif_match(patients_data, motif):
    ms_count = sum([has_motif_match(e, motif) for e in patients_data['MS']])
    hc_count = sum([has_motif_match(e, motif) for e in patients_data['HC']])
    return {'score': ms_count - hc_count,
            'ms_count': ms_count,
            'hc_count': hc_count}

def has_motif_match(patient_seq_info, motif):

    sliced_indices = []
    sliced_motif_array = []
    # Make sure the explict motif subsequences are sorted by smallest to largest positions relative
    # to the cdr
    sorted_motif = sorted(motif, key=lambda x: x[1])
    for (explicit_motif_subseq, rel_position) in sorted_motif:

        # Negative `rel_position' indicates that the explicit motif subsequence occurs to
        # the left of the cdr, and positive indicates that starting point is to the right. In
        # either case, adding in rel_position to 'max_cdr_start' will ensure that we obtain the
        # correct index.
        index = patient_seq_info['max_cdr_start'] + rel_position
        sliced_indices.extend(range(index, index + len(explicit_motif_subseq)))
        sliced_motif_array.extend([AMINO_ACID_CODES.get(aa, len(AMINO_ACID_CODES) + 1)
                                   for aa in explicit_motif_subseq])
    return np.any(np.all(patient_seq_info['seq_matrix'][:, sliced_indices] == sliced_motif_array,
                         axis=1))



def load_patients_data(patient_dir='example_patient_files', delimiter='_', diagnosis_index=1,
                       patient_id_index=2, ms_marker='MS', hc_marker='HC'):

    patients_data = {ms_marker: [], hc_marker: []}
    for name in os.listdir(patient_dir):
        name_list = name.replace('.fa', '').split(delimiter)
        diagnosis = name_list[diagnosis_index]
        if diagnosis == ms_marker:
            diagnosis = 'MS'
        elif diagnosis == hc_marker:
            diagnosis = 'HC'
        else:
            continue

        patient_id = name_list[patient_id_index]
        fname = os.path.join(patient_dir, name)
        patients_data[diagnosis].append(load_patient_sequences(patient_id, fname))
    return patients_data

def load_patient_sequences(patient_id, fname, delimiter=';', cdr_index=-1):

    num_seqs = 0
    max_cdr_start = 0
    max_pos_cdr_len = 0

    with open(fname) as f:
        for line in f:
            if line.startswith('>'):
                cdr = line.split(delimiter)[cdr_index].strip()
                continue
            else:
                if cdr not in line:
                    continue

                num_seqs += 1
                max_cdr_start = max(max_cdr_start, line.find(cdr))
                max_pos_cdr_len = max(max_pos_cdr_len, len(line.strip()) - line.find(cdr))

    patient_matrix = np.zeros((num_seqs, max_cdr_start + max_pos_cdr_len), dtype=np.uint8)

    seq_index = 0
    with open(fname) as f:
        for line in f:
            if line.startswith('>'):
                cdr = line.split(delimiter)[cdr_index].strip()
                continue
            else:
                if cdr not in line:
                    continue

            sequence = line.strip()
            cdr_start = line.find(cdr)
            for i, aa in enumerate(sequence[:cdr_start]):
                # Ensure the the index is position relative both cdr_start and max_cdr_start for
                # all sequences in patient
                index = max_cdr_start + (i - cdr_start)
                amino_acid_val = AMINO_ACID_CODES.get(aa, len(AMINO_ACID_CODES) + 1)
                patient_matrix[seq_index][index] = amino_acid_val
            for i, aa in enumerate(sequence[cdr_start:]):
                index = max_cdr_start + i
                amino_acid_val = AMINO_ACID_CODES.get(aa, len(AMINO_ACID_CODES) + 1)
                patient_matrix[seq_index][index] = amino_acid_val

            seq_index += 1
    patient_sequence_info = {'patient_id': patient_id,
                             'seq_matrix': patient_matrix,
                             'max_cdr_start': max_cdr_start,
                             'max_pos_cdr_len': max_pos_cdr_len}

    return patient_sequence_info


if __name__ == '__main__':
    total_doms = sorted(run_haystack(), key=lambda x:x[1])
    for dom, count in total_doms:
        print "DOM %s occurs in %d patients"%(str(dom), count)


