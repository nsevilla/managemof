### Short step-by-step to create MOF files ready for upload (for Y3A1). This has to be run within the DESDM system to have access to /archive_data

1. Create a table with unitname, reqnum, attnum for the Y3 objects in the database for each tile. From easyaccess for instance:

select unitname, reqnum, attnum from proctag t, pfw_attempt p where t.tag = 'Y3A1_MOF' and p.id = t.pfw_attempt_id; > FILE_WITH_TABLE.tab (you may want to remove the header of this file for subsequent steps, I did it by hand using sed -i 1d FILE_WITH_TABLE.tab)

2. Run shell script link_to_moftiles.sh FILE_WITH_TABLE.tab. This will create a link for each combination of  unitname, reqnum, attnum in the file created in step 1. Directory path has to be modified in the script (e.g. /archive_data/desarchive/OPS/multiepoch/Y3A2_MOF). All links will be at a directory called LINKDIR for the purposes of this README.

3. Run python script download_ids_for_tiles.py -d LINKDIR -i IDSDIR. pyfits 3.3+8 and despydb 2.0.2+1 are necessary (at least, these versions work). This will download a fits file per tile present in the LINKDIR directory containing the COADD_OBJECT_ID of only those objects present in the Y3A2_COADD_OBJECT_SUMMARY table (this is hardcoded for now). Files are downloaded to IDSDIR.

4. Run python script flatten_match_mof.py -d LINKDIR -i IDSDIR -o UPLOADDIR. For each file in IDSDIR, create a file in UPLOADDIR for upload, containing the flattened data in the files from the LINKDIR directory (the MOF files) but only for those objects included in the corresponding ID file in IDSDIR.  

5. Merge the output files into large files to optimize upload. E.g. using ftmerge @filelist.txt upload.fits

6. Ingest into the DB. E.g. easyaccess -c "load_table upload.fits"
