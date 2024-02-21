import scanpy as sc
import muon

def load_h5mu(h5mu_file):
    return muon.read_h5mu(h5mu_file)


def compute_umap(adata):
    '''
    '''
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
    sc.tl.umap(adata)

def compute_tsne(adata):
    '''
    '''