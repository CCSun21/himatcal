"""utils for covalent"""
import covalent as ct
from ase.io import write


def extract_result(dispatch_id, return_atoms=True):    # sourcery skip: avoid-builtin-shadow
    """
    extract the result of a covalent job after Gaussian calculation
    """
    import os
    import pickle as pkl

    from monty.os import cd
    results = ct.get_result(dispatch_id)
    dispatch_id = results.dispatch_id
    # create a folder with the dispatch_id
    os.makedirs(dispatch_id, exist_ok=True)
    with cd(dispatch_id):
        with open(f'{dispatch_id}.pkl','wb') as f:
            pkl.dump(results,f)
        try:
            if return_atoms:
                return extract_gaussian_result(results, dispatch_id)
            else:
                extract_gaussian_result(results, dispatch_id)
                return results # return covalent results
        except Exception:
            if results.status.__eq__(_Status__value='RUNNING'):
                print('Job is still running')
            with open('log.yaml','w') as f:
                f.write(f'dispatch_id: {dispatch_id}\n')
                f.write(f'job_status: {str(results.status.__str__)}\n')


def extract_gaussian_result(results, id):
    final_scf_energy = results.result['attributes']['final_scf_energy']/27.211385 # convert the uint from a.u. to eV
    parameters = results.result['parameters']
    input_atoms = results.result['trajectory'][0]
    final_atoms = results.result['atoms']
    write('input.xyz',input_atoms)
    write('final.xyz',final_atoms)
    with open('log.yaml','w') as f:
        f.write(f'dispatch_id: {id}\n')
        f.write(f'status: \'{str(results.status.__str__)}\'\n')
        f.write(f'final_scf_energy (a.u.): {final_scf_energy}\n')
        f.write(f'parameters: {parameters}\n')
    return final_atoms


def get_result(dispatch_id, job_dir, wait=False):
    """
    extract result in the job_dir
    """
    from monty.os import cd
    results = ct.get_result(
        dispatch_id,
        wait=wait
    )
    if results.status.__eq__("COMPLETED") is True:
        with cd(job_dir):
            extract_result(dispatch_id)
            final_scf_energy = (
                results.result["attributes"]["final_scf_energy"] / 27.211385
            )  # convert to a.u.
            atoms = results.result["atoms"]
            write("final.xyz", atoms)
            print(f"final_scf_energy (a.u.): {final_scf_energy}\n")
        return atoms
    else:
        print("The calculation is not completed.")
        return None