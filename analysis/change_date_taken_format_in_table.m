% change dateTaken format and move variables around
brx_num = 2;
tablename = strcat(projectdir,'\',num2str(brx_num),'\',experimentnumber,"_",num2str(brx_num),".mat");
load(tablename)
dtable.dateTaken.Format = 'ddMMyyyy HHmmss';
dtable = movevars(dtable,'markerArea','After','imageName');
dtable = movevars(dtable,'pressure','Before','center1');
dtable = movevars(dtable,'flow','Before','center1');

save(tablename,'dtable');
