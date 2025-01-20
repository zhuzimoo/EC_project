from mebocost import mebocost
import scanpy as sc

tissue_name = "intestine"
adata = sc.read_h5ad("./raw_10000/intestine_10000.h5ad")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


mebo_obj = mebocost.create_obj(
                        adata = adata,
                        group_col = ['ct'],
                        met_est = 'mebocost',
                        config_path = '/temp_work/ch245921/zimozhu/MEBOCOST/mebocost.conf',
                        exp_mat=None,
                        cell_ann=None,
                        species='human',
                        met_pred=None,
                        met_enzyme=None,
                        met_sensor=None,
                        met_ann=None,
                        scFEA_ann=None,
                        compass_met_ann=None,
                        compass_rxn_ann=None,
                        gene_network=None,
                        gmt_path=None,
                        cutoff_exp='auto', ## automated cutoff to exclude lowly ranked 25% sensors across all cells
                        cutoff_met='auto', ## automated cutoff to exclude lowly ranked 25% metabolites across all cells
                        cutoff_prop=0.25, ## at lease 25% of cells should be expressed the sensor or present the metabolite in the cell group (specified by group_col)
                        sensor_type=['Receptor', 'Transporter', 'Nuclear Receptor'],
                        thread=8
                        )
print("mebocost created")

commu_res = mebo_obj.infer_commu(
                                n_shuffle=1000,
                                seed=12345,
                                Return=True,
                                thread=None,
                                save_permuation=False,
                                min_cell_number = 1
                            )

mebocost.save_obj(obj = mebo_obj, path = f'./comb_result/mebocost/pk_file/{tissue_name}_commu.pk')
commu_res = mebo_obj.commu_res.copy()
commu_res = commu_res[commu_res['permutation_test_fdr']<=0.05]
commu_res.to_csv(f'./comb_result/mebocost/csv_file/{tissue_name}_communication_result.tsv', sep = '\t', index = None)
