from gcloud import datastore
import requests
import time
ds = datastore.Client(project='andersen-lab')

url = "https://docs.google.com/spreadsheets/d/1V6YHzblaDph01sFDI8YK_fP0H7sVebHQTXypGdiQIjI/pub?gid=0&single=true&output=tsv"

gs = requests.get(url).text.encode("utf-8").splitlines()

gs = [str(x, 'utf-8').strip().split("\t") for x in gs]
gs = [x for x in gs if x[2]]

gs = [(x[0], x[2], x[3]) for x in gs]
WI_ISOTYPE = {}
WI_STRAIN = {}
for strain, isotype, prev_names in gs:
    if prev_names != "NA":
        prev_names = prev_names.split("|")
        for p in prev_names:
            if p:
                WI_ISOTYPE[p] = isotype
                WI_STRAIN[p] = strain
    if strain and isotype:
        WI_ISOTYPE[strain] = isotype
        WI_STRAIN[strain] = strain

    if isotype:
        WI_ISOTYPE[isotype] = isotype

exclude_indices = ['most_abundant_sequence',
                   'fastqc_per_base_sequence_quality_data',
                   'fastqc_per_tile_sequence_quality_data',
                   'fastqc_per_sequence_quality_scores_data',
                   'fastqc_per_base_sequence_content_data',
                   'fastqc_per_sequence_gc_content_data',
                   'fastqc_per_base_n_content_data',
                   'fastqc_sequence_length_distribution_data',
                   'fastqc_sequence_duplication_levels_data',
                   'fastqc_overrepresented_sequences_data',
                   'fastqc_adapter_content_data',
                   'fastqc_kmer_content_data',
                   'fastqc_error']


def query_item(kind, filters=None, projection=()):
    # filters:
    # [("var_name", "=", 1)]
    query = ds.query(kind=kind, projection = projection)
    if filters:
        for var, op, val in filters:
            query.add_filter(var, op, val)
    return query.fetch()

def get_item(kind, name):
    return ds.get(ds.key(kind, name))

def update_item(kind, name, **kwargs):
    item = get_item(kind, name)
    if item is None:
        m = datastore.Entity(key=ds.key(kind, name),
                             exclude_from_indexes=exclude_indices)
    else:
        m = datastore.Entity(key=ds.key(kind, name),
                             exclude_from_indexes=exclude_indices)
        m.update(dict(item))
    for key, value in kwargs.items():
        if type(value) == str:
            m[key] = value
        elif type(value) == list:
            if key in m:
                m[key] += value
            else:
                m[key] = value
            m[key] = list(set(m[key]))
        # If date created of file is earlier
        elif key == 'date_created' and item:
            vtimestamp = time.mktime(value.timetuple())
            dstimestamp = time.mktime(m['date_created'].timetuple())
            if vtimestamp < dstimestamp:
                m[key] = value

        else:
            m[key] = value
    if 'fq_profile_count' in m:
        m['fq_profile_count'] += 1
    else:
        m['fq_profile_count'] = 1
    ds.put(m)




fastqs = query_item('fastq', filters = [['strain_type', '=', 'WI'], ['use', '=', True]])

for fq in fastqs:
    if 'original_strain' in fq.keys():
        if fq['original_strain'] in WI_STRAIN.keys():
            fq['strain'] = WI_STRAIN[fq['original_strain']]
        if fq['original_strain'] in WI_ISOTYPE.keys():
            fq['isotype'] = WI_ISOTYPE[fq['original_strain']]
            print([fq.key.name, fq['isotype'], fq['strain'], fq['original_strain']])
    if 'seq_folder' in fq.keys():
        if fq['seq_folder'] != "original_wi_seq":
            if fq['library'] != fq['barcode'].replace("+", ""):
                print(fq['library'] + "=>" + fq['barcode'].replace("+", ""))
                fq['library'] = fq['barcode'].replace("+", "")
    update_item('fastq', fq.key.name, **fq)
