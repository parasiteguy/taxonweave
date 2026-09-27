from taxonweave.worms import (
    resolve_name,
    get_accepted_record,
    get_classification,
    flatten_classification,
    get_synonyms
)


# Try an unaccepted historical name
species = "Mercierella enigmatica"


# --------------------------------------------------
# RESOLVE USER-SUPPLIED NAME
# --------------------------------------------------

record = resolve_name(species)


if record is None:

    print("No WoRMS record found.")

else:

    print("\nQUERY\n")

    print("Name entered:", species)
    print("WoRMS match:", record.get("scientificname"))
    print("Status:", record.get("status"))
    print("AphiaID:", record.get("AphiaID"))


    # --------------------------------------------------
    # RESOLVE TO ACCEPTED TAXON
    # --------------------------------------------------

    accepted = get_accepted_record(record)

    accepted_aphia_id = accepted.get("AphiaID")

    print("\nACCEPTED TAXON\n")

    print("Scientific name:", accepted.get("scientificname"))
    print("Authority:", accepted.get("authority"))
    print("Status:", accepted.get("status"))
    print("Rank:", accepted.get("rank"))
    print("AphiaID:", accepted_aphia_id)


    # --------------------------------------------------
    # CLASSIFICATION
    # --------------------------------------------------

    classification = get_classification(accepted_aphia_id)

    flattened = flatten_classification(classification)

    print("\nCLASSIFICATION\n")

    for taxon in flattened:

        print(
            f"{taxon['rank']}: "
            f"{taxon['scientificname']} "
            f"(AphiaID {taxon['AphiaID']})"
        )


    # --------------------------------------------------
    # SYNONYMS
    # --------------------------------------------------

    synonyms = get_synonyms(accepted_aphia_id)

    print("\nSYNONYMS\n")

    if not synonyms:

        print("No synonyms found.")

    else:

        for synonym in synonyms:

            print(
                f"{synonym.get('scientificname')} "
                f"{synonym.get('authority')} "
                f"[AphiaID {synonym.get('AphiaID')}]"
            )