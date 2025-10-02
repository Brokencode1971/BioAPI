from flask import Flask, render_template, request, jsonify
import requests

app = Flask(__name__)

def safe_get(d, key, default=None):
    return d.get(key, default) if isinstance(d, dict) else default

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/search", methods=["POST"])
def search():
    protein_id = request.form.get("query")
    if not protein_id:
        return jsonify({"error": "No ID provided."}), 400

    # Fetch UniProt JSON
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.json"
    response = requests.get(url)
    if response.status_code != 200:
        return jsonify({"error": "Protein not found."}), 404

    data = response.json()
    app.logger.info(f"UniProt API response for {protein_id}: {data}")

    try:
        # ——— OVERVIEW EXTRACTION ———
        overview = {}
        overview['accession'] = data.get('primaryAccession')
        overview['recommended_name'] = safe_get(
            safe_get(
                safe_get(data, 'proteinDescription', {}),
                'recommendedName',
                {}
            ),
            'fullName',
            {}
        ).get('value', 'N/A')
        overview['alternative_names'] = [
            safe_get(alt, 'fullName', {}).get('value', 'N/A')
            for alt in safe_get(data.get('proteinDescription', {}), 'alternativeNames', [])
            if isinstance(safe_get(alt, 'fullName'), dict)
        ]

        # Gene names + synonyms
        gene_names = []
        for gene in data.get('genes', []):
            if safe_get(gene, 'geneName'):
                gene_names.append(safe_get(gene['geneName'], 'value', 'N/A'))
            for syn in gene.get('synonyms', []):
                gene_names.append(safe_get(syn, 'value', 'N/A'))
        overview['gene_names'] = gene_names

        overview['organism'] = safe_get(data.get('organism', {}), 'scientificName', 'N/A')
        overview['protein_existence'] = safe_get(data.get('proteinExistence', {}), 'value')
        overview['entry_version'] = data.get('entryVersion')
        overview['sequence_version'] = safe_get(data.get('sequence', {}), 'version')

        # Function summary
        function_texts = []
        for comment in data.get('comments', []):
            if comment.get('commentType') == 'FUNCTION':
                for txt in comment.get('texts', []):
                    function_texts.append(safe_get(txt, 'value', 'N/A'))
        overview['function'] = function_texts

        # ——— Sequence & Basic Panels ———
        sequence = safe_get(data.get('sequence', {}), 'value', 'N/A')

        pdb_ids = [xref['id'] for xref in data.get('uniProtKBCrossReferences', []) if xref.get('database') == 'PDB']
        pdb_text = None
        if pdb_ids:
            pdb_url = f"https://files.rcsb.org/download/{pdb_ids[0]}.pdb"
            resp = requests.get(pdb_url)
            if resp.status_code == 200:
                pdb_text = resp.text

        # Cross-links
        url_templates = {
            'Ensembl':  'https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={id}',
            'KEGG':     'https://www.genome.jp/dbget-bin/www_bget?hsa:{id}',
            'Reactome': 'https://reactome.org/PathwayBrowser/#/{id}',
            'GeneID':   'https://www.ncbi.nlm.nih.gov/gene/{id}',
            'Pfam':     'https://pfam.xfam.org/protein/{id}',
            'InterPro': 'https://www.ebi.ac.uk/interpro/entry/UniProt/{protein_id}'
        }
        cross_links = []
        for xref in data.get('uniProtKBCrossReferences', []):
            db = xref.get('database')
            _id = xref.get('id')
            if db in url_templates:
                cross_links.append({
                    'db': db,
                    'id': _id,
                    'url': url_templates[db].format(id=_id, protein_id=protein_id)
                })

        # ——— New Detailed Panels ———
        # Activity
        activity = {'catalytic_activity': [], 'regulation': []}
        for comment in data.get('comments', []):
            ct = comment.get('commentType')
            if ct == 'CATALYTIC ACTIVITY':
                rxn = comment.get('reaction', {})
                activity['catalytic_activity'].append(rxn.get('name'))
            elif ct == 'ACTIVITY REGULATION':
                for txt in comment.get('texts', []):
                    activity['regulation'].append(safe_get(txt, 'value', ''))

        # Localization
        localization = []
        for comment in data.get('comments', []):
            if comment.get('commentType') in ['SUBCELLULAR LOCATION', 'TISSUE SPECIFICITY']:
                if 'subcellularLocations' in comment:
                    for loc in comment['subcellularLocations']:
                        localization.append(loc['location']['value'])
                else:
                    for txt in comment.get('texts', []):
                        localization.append(safe_get(txt, 'value', ''))

        # Modifications & Domains & Sites
        modifications = []
        domains_sites = []
        for feat in data.get('features', []):
            t = feat.get('type')
            if t in ['Glycosylation', 'Disulfide bond', 'PTM']:
                modifications.append(f"{t} at {feat['location']['start']['value']}")
            elif t in ['Domain', 'Active site', 'Signal', 'Propeptide', 'Chain']:
                domains_sites.append(f"{t}: {feat.get('description', '')} ({feat['location']['start']['value']}-{feat['location']['end']['value']})")

        # Keywords
        keywords = [kw['name'] for kw in data.get('keywords', [])]

        # References
        references = []
        for ref in data.get('references', []):
            cit = ref.get('citation', {})
            title = cit.get('title')
            pm = next((cr['id'] for cr in cit.get('citationCrossReferences', []) if cr['database']=='PubMed'), None)
            references.append(f"{title} (PubMed:{pm})")

        # ——— Build Response ———
        return jsonify({
            'overview': overview,
            'sequence': sequence,
            'structure': {'pdb_ids': pdb_ids, 'pdb_text': pdb_text},
            'links': cross_links,
            'activity': activity,
            'localization': localization,
            'modifications': modifications,
            'domains_sites': domains_sites,
            'keywords': keywords,
            'references': references
        })

    except Exception as e:
        app.logger.error(f"Error in /search: {e}")
        return jsonify({"error": f"Server error: {e}"}), 500

@app.route("/autocomplete")
def autocomplete():
    query = request.args.get("q") or ""
    url = (
        "https://rest.uniprot.org/uniprotkb/search"
        f"?query={query}"
        "&fields=accession,protein_name"
        "&format=json"
        "&size=5"
    )
    resp = requests.get(url)
    data = resp.json() if resp.status_code == 200 else {}
    results = []
    for item in data.get("results", []):
        results.append({
            'id': item.get('primaryAccession'),
            'name': item.get('proteinDescription', {})
                       .get('recommendedName', {})
                       .get('fullName', {})
                       .get('value', '')
        })
    return jsonify(results)

if __name__ == "__main__":
    app.run(debug=True)
