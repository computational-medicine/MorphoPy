function [filename_unpacked2] = unzip_temp (filename)
% [tmpname] = unzip_temp (filename);
%
% Unpacks a zipped file in .gz format to a temporary location. The original file
% is preserved.
%
%  filename       Name of the original file that is to be unzipped
%
%  tmpname       Name of the temporary file that is generated. This file can
%                 (and should!) safely be deleted after use.
%
% Example:
%  t = unzip_temp ('~/foo.nii.gz');
%  % ... process t ...
%  delete (t);

   % quote any spaces in the filename (to avoid using quotes in the commands)
   %escaped = strrep (filename, ' ', '\ ');
   %tempdir
   %pause
   %filename
   %pause
   k10=strfind(filename, '\');
   if isempty(k10)
       filename_exact=filename;
   else
       filename_exact=filename(max(k10)+1:end);
   end
   %pause
   filename_unpacked=strrep(filename, '.gz', '');
   %filename_unpacked
   %pause
   escaped=strrep (filename, '\', '\\');
   %escaped=filename;
   escaped = ['"' escaped '"'];
   %pause
   % first shot at parsing the name; it should end with .gz to indicate that
   % this is indeed a zipped file
   [a, name, ext] = fileparts (filename);
   if ~strcmp (ext, '.gz'),
      error ('File "%s" does not have .gz extension', filename);
   end;
   
   % parse the name *again*, to recover the original extension
   [a, stem, real_ext] = fileparts (name);
   
   % don't use spaces in stem, rather replace with underscores
   stem = strrep (stem, ' ', '_');
   
   % query the environment for clues as to where temporary files should be
   % stored
   sys_tmp = getenv ('TEMP');
   if isempty (sys_tmp),
         sys_tmp = '/tmp';
   end;
   
   % hash in the user name, to avoid cluttering the namespace
   usr_name = getenv ('USER');   
   tmp_dir = fullfile (sys_tmp, usr_name);
   
   % make sure that this directory exists
   if ~exist (tmp_dir, 'dir'),
      mkdir (tmp_dir);
   end;
   
%    % generate a temporary name; unfortunately the mktemp command does not (yet)
%    % support suffices
%    [err_code, tmp_stem] = system (sprintf ('mktemp -p %s %s_XXXXXX', ...
%                                                             tmp_dir, stem));
%                                                         
%    if err_code,
%       error ('mktemp returned %d for directory "%s" and name "%s"', ...
%                                                    err_code, tmp_dir, stem);
%    end;
%    
%    % remove newline from the output (!) and quote spaces here too
%    tmp_stem = strrep (strrep (tmp_stem, sprintf ('\r'), ''), ...
%                                                             sprintf ('\n'), '');
%    
%    % compose the final name of the file, inserting back the real extension
%    tmpname = horzcat (tmp_stem, real_ext);
    
    % temporary filename
%     tmpname = [tempname '.nii'];
%     tmpname_out = tmpname;
%     tmpname = ['"' tmpname '"'];
%     tmpname
%     pause
%     tmpname=strrep (tmpname, '\', '\\');
    
    tmpname2 = [tempdir filename_exact];
    filename_unpacked2=strrep(tmpname2, '.gz', '');
    %tmpname_out2 = tmpname2;
    %tmpname2 = ['"' tmpname2 '"'];
    %tmpname2
    %pause
    %tmpname2=strrep (tmpname2, '\', '\\');

%    % rename temporary file (which already exists) to final name
%    [err_code] = system (sprintf ('mv %s %s', tmp_stem, tmpname));
%    if err_code,
%       error ('mv returned %d for source "%s" and target "%s"', ...
%                                                 err_code, tmp_stem, tmpname);
%    end;
   
   % decompress the file into the temporary location
    %    v = sprintf ('gunzip %s > %s', escaped, tmpname)
%     if ispc
%         [err_code] = system(['7za x ' tmpname ' > ']);
%     else
%         [err_code] = system (sprintf ('zcat %s > %s', escaped, tmpname));
%     end;
%     tempname
%     pause
%     1
%     sprintf (escaped, tmpname)
%     pause
%     2
%    sprintf ('C:\\Program Files (x86)\\gzip\\gzip')
%    pause
%    3
%    sprintf ('C:\\Program Files (x86)\\gzip\\gzip -cd %s > %s', escaped, tmpname)
%    pause
%    4
%     sprintf ('%s %s',escaped, tmpname)
%    pause
%    %[err_code] = system (sprintf ('C:\\Program Files (x86)\\gzip\\gzip -cd %s > %s', escaped, tmpname));
%    5
%    %['C:\Program Files (x86)\gzip\gzip -cd ' escaped ' ' tmpname]
%    ['C:\Program Files (x86)\7-zip\7z e ' escaped ' ' tmpname]
%    pause
%(sprintf('"C:\\Program Files (x86)\\7-zip\\7z" e %s %s', escaped, tmpname))
%pause
%((['"C:\\Program Files (x86)\\7-zip\7z" x ' escaped ' -o' tempdir]))
%pause

%tmpname
%pause
   [err_code] = system ((['"C:\\Program Files (x86)\\7-zip\\7z" x ' escaped ' -o' tempdir]));
   %escaped
   %[err_code] = system (sprintf ('C:\\Program Files (x86)\\gzip\\gzip -cd %s > %s', escaped, tempdir));
   %[err_code] = system (sprintf ('C:\\Program Files (x86)\\gzip\\gzip'));
   if err_code,
      error ('zcat returned %d for source "%s" and target "%s"', ...
                                                err_code, filename, tmpname);
                                            %error ('zcat returned %d for source "%s"', ...
                                            %    err_code, filename);
   end;
%10

end