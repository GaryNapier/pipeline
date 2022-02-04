# tb_data class

# Read in data from Jody's json results files and store as a simplified version

# from tqdm import tqdm
import json
import pathogenprofiler as pp
from tqdm import tqdm


class tb_data:
    pass

    # WHY NO 'self' HERE?
    # myVar = 10

    # Example usage:
    # filename = '/mnt/storage7/jody/tb_ena/tbprofiler/freebayes/results/'
    # genes = ('ahpC', 'katG', 'fabG1')
    # all_data = tb_data(filename, suffix, meta_dict, genes).all_data

    # def __init__(self, filename, suffix, meta_dict, genes):

    #     standardise_drtype = {
    #     "Sensitive":"Sensitive",
    #     "Pre-MDR":"Pre-MDR-TB",
    #     "HR-TB":"Pre-MDR-TB",
    #     "RR-TB":"Pre-MDR-TB",
    #     "MDR":"MDR-TB",
    #     "MDR-TB":"MDR-TB",
    #     "Pre-XDR":"Pre-XDR-TB",
    #     "Pre-XDR-TB":"Pre-XDR-TB",
    #     "XDR":"XDR-TB",
    #     "XDR-TB":"XDR-TB",
    #     "Other":"Other"
    #     }

    #     # Read in all the json data for samples with given genes only (>0.7 freq); combine with the metadata for those samples
    #     self.all_data = {}
    #     for samp in tqdm(meta_dict):
    #         self.tmp = []
    #         # self.tmp = []
    #         try:
    #             # Open the json file for the sample, skip if can't find
    #             tmp_data = json.load(open(pp.filecheck("%s/%s%s" % (filename, samp, suffix))))
    #             # self.tmp_data = json.load(open(pp.filecheck("%s/%s%s" % (file, samp, suffix))))
    #             # print("LOADED SAMP", samp)
    #         except:
    #             pass

    #         for var in tmp_data["dr_variants"] + tmp_data["other_variants"]:
    #         # for var in self.tmp_data["dr_variants"] + self.tmp_data["other_variants"]:
    #             if var['gene'] not in genes: continue
    #             if var['freq'] < 0.7: continue
    #             tmp.append(var)
    #             # self.tmp.append(var)
    #         # If sample meets criteria above, then tmp list will not be empty
    #         if len(tmp) > 0:
    #         # if len(self.tmp) > 0:
    #             # Create empty dict for next two entries
    #             self.all_data[samp] = {}
    #             # Append metadata dict
    #             self.all_data[samp]['metadata'] = {
    #                 'wgs_id':samp,
    #                 'inh_dst':meta_dict[samp]['isoniazid'],
    #                 'main_lin': tmp_data['main_lin'],
    #                 # 'main_lin': self.tmp_data['main_lin'],
    #                 'sublin':tmp_data['sublin'],
    #                 # 'sublin':self.tmp_data['sublin'],
    #                 'country_code':meta_dict[samp]['country_code'],
    #                 'drtype':standardise_drtype[tmp_data['drtype']]
    #                 # 'drtype':standardise_drtype[self.tmp_data['drtype']]
    #                 }
    #             # Append mutations dict
    #             self.all_data[samp]['mutations'] = tmp

    #     # Find mixed samples
    #     mixed_samp_lin_dict = {}
    #     for samp in self.all_data:
    #         lin = self.all_data[samp]['metadata']['main_lin']
    #         sublin = self.all_data[samp]['metadata']['sublin']
    #         if ";" in lin+sublin:
    #             mixed_samp_lin_dict[samp] = {'lin': lin, 'sublin': sublin}

    #     # Remove from data
    #     for samp in list(self.all_data):
    #         if samp in mixed_samp_lin_dict:
    #             del self.all_data[samp]














