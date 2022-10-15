
 function [numpolys,polys] = fixed_upright(n)
% Generate and enumerate all the fixed polyiomonds of size n triangles, startng with an inverted iamond.
% Returns the total number of fixed polyiomonds (numpolys) and polys, a
% cell array of size numpolys, where each cell is a nx2 array of (x,y)
% positions of each triangle in the polyiomond.

% Uses algorithm from https://www.sciencedirect.com/science/article/pii/0012365X81902375?via%3Dihub
% Inspired by MATLAB code for polyominoes given by Aaron T. Becker and Yitong Lu (https://in.mathworks.com/matlabcentral/fileexchange/75698-polyenumdraw?s_tid=prof_contriblnk)

polyList = zeros(n,1); % which triangles in adjtriangles are used
adjtriangles = zeros(2*n,2);

polyList(1) = 1; % set first triangle was used

adjtriangles(1:2,:) = [0,0;1,0];

numtriangles = 1;

% The simplest implementation involves adding one triangle at a time.
% Beginning with an initial triangle, number the adjacent triangles clockwise from the top, 1, 2, 3, and 4.
% Do not number any triangle that is on a lower row, or left of the triangle on the same row.
% This is the version described by Redelmeier.

numAdj = 2;
[numpolys,polys] = recursPolyBuild(adjtriangles, polyList, numtriangles,numAdj,n);

    function [numpolys,polys] = recursPolyBuild(adjtriangles, polyList, numtriangles, numAdj,n)
        numpolys = 0;
        polys ={};
        if numtriangles+1 == n
            numpolys = numAdj-polyList(numtriangles);
            for i = polyList(numtriangles)+1 : numAdj
                
                polys{end+1} = adjtriangles([polyList(1:n-1);i],:); %#ok<AGROW>
               
            end
            return;
        end
        
        % pick a number between polyList(numtriangles)+1 and numAdj,
        for i = polyList(numtriangles)+1 : numAdj
            adjtrianglesi = adjtriangles;
            polyListi = polyList;
            numAdji = numAdj;
            b = adjtriangles(i,:);
            % add a triangle at that location.
            polyListi(numtriangles+1) = i;
            
            % Number the unnumbered adjacent triangles, starting with 5, remove ones already in adjList
           if mod(b(1,1)+b(1,2),2)~=0  %If sum of coordinates is odd,triangle is inverted and will have a restriction in adjacency below.
            a = b + [0,1];
            if ~any(adjtriangles(:,1)==a(1) & adjtriangles(:,2)==a(2))
                numAdji = 1 + numAdji;
                adjtrianglesi(numAdji,:) = a;
            end
           else
               a = b + [0,-1];
            if a(2)>0 || (a(2)==0 && a(1)>0)  % don't add triangles below or to the left of starting triangle @ (0,0)
                if ~any(adjtriangles(:,1)==a(1) & adjtriangles(:,2)==a(2))
                    numAdji = 1 + numAdji;
                    adjtrianglesi(numAdji,:) = a;
                end
            end           
           end
            a = b + [1,0];
            if ~any(adjtriangles(:,1)==a(1) & adjtriangles(:,2)==a(2))
                numAdji = 1 + numAdji;
                adjtrianglesi(numAdji,:) = a;
            end
            
            
            
            
            a = b + [-1,0];
            if a(2)>0 || (a(2)==0 && a(1)>0)  % don't add triangles below or to the left of starting triangle @ (0,0)
                if ~any(adjtriangles(:,1)==a(1) & adjtriangles(:,2)==a(2))
                    numAdji = 1 + numAdji;
                    adjtrianglesi(numAdji,:) = a;
                end
            end
            [numpolysN,polysN] = recursPolyBuild(adjtrianglesi, polyListi, numtriangles+1, numAdji,n);
          numpolys = numpolys +numpolysN;
             polys = [polys,polysN]; %#ok<AGROW>
            
            % Then, pick a number larger than the previously picked number,
            % and add that triangle
            % Continue picking a number larger than the number of the
            % current triangle, adding that triangle, and then numbering the
            % new adjacent triangles. When n triangles have been created, an
            % n-iomond has been created.
        end
       
    end
    
 end 
  
