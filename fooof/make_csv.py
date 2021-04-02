def make_csv(csv_path, extension):
    """Function to insert the number of epochs to include in analysis into csv.
    Number of epochs is calculated by comparing the number of epochs available
    for each subject and including the minimum amount.

    Parameters
    ----------
    csv_path : str,
        path to the csv containing information on subjects to include

    extension : str,
        file extension of meg files (e.g. .asc)

    Returns
    -------
    None
    saves the extended csv to the same directory where found old csv
    (i.e. overwrites old csv)

    """

    df = pd.read_csv(csv_path, delimiter =  ',')

    nr_epochs = []
    for index, row in df.iterrows():
        asc_paths = find_paths(main_dir=row['Path'],
                            subject=row['Case_ID'],
                            extension=extension,
                            timepoint=row['MM'],
                            atlas=row['Atlas'])
        #store nr of epochs available for each subject
        nr_epochs.append(len(asc_paths))

    #find smallest number of epochs available
    min_nr_epochs = min(nr_epochs)

    #add start and stop epochs to df
    df['Start'] = np.repeat(1,len(df['Path']))
    df['End'] = np.repeat(min_nr_epochs, len(df['Path']))

    #save new csv file that includes the epochs to analyse
    df.to_csv(csv_path)