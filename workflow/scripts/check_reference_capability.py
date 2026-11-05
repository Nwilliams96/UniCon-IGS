"""
Check whether the installed RESCRIPt environment supports the requested reference release.
This preflight check runs before read processing. It inspects capabilities only; reference acquisition and preparation remain separate Snakemake rules.
"""

import json  #Provides JSON serialization for the capability record.
from pathlib import Path  #Provides platform-independent output-path handling.


#This function validates the requested database release and records the compatible software versions.
def check_capability(database, version, output):
    import qiime2  #Provides the installed QIIME 2 version.
    import rescript  #Provides the installed RESCRIPt version and plugin.
    import sklearn  #Provides the scikit-learn version associated with the classifier environment.
    from qiime2.sdk import PluginManager  #Provides access to registered QIIME 2 plugin actions.

    action_name = "get_silva_data" if database == "SILVA" else "get_pr2_data"  #Selects the RESCRIPt acquisition action for the requested database.
    action = PluginManager().plugins["rescript"].actions[action_name]  #Retrieves the selected action from the installed RESCRIPt plugin.
    supported_type = action.signature.parameters["version"].qiime_type  #Retrieves the allowed values for the action’s version parameter.
    if version not in supported_type:  #Stops when the requested release is not supported by this environment.
        raise ValueError(
            f"Installed RESCRIPt {rescript.__version__} cannot download {database} {version}. "  #Reports the incompatible database release.
            f"Supported versions: {supported_type}. Select a supported training release, "  #Lists the releases accepted by the installed action.
            f"supply {database}_download_link for a compatible pre-trained classifier, "  #Suggests using a compatible pre-trained classifier.
            "or use an environment with support for the requested release."  #Suggests switching to a compatible environment.
        )
    path = Path(output)  #Converts the declared output into a Path object.
    path.parent.mkdir(parents=True, exist_ok=True)  #Creates the output directory when it does not already exist.
    path.write_text(json.dumps(dict(database=database, version=version,  #Records the validated database and release.
                                   qiime2=qiime2.__version__, rescript=rescript.__version__,  #Records the QIIME 2 and RESCRIPt versions.
                                   sklearn=sklearn.__version__), indent=2) + "\n")  #Records scikit-learn and writes formatted JSON with a final newline.


#This block runs the capability check when the script is executed by Snakemake.
if "snakemake" in globals():  #Avoids automatic execution when the module is imported independently.
    check_capability(snakemake.params.database, snakemake.params.version, snakemake.output[0])  #Checks the requested release and writes the capability marker.
