classdef circularGraph_noButtons_colorGradient < handle
    % CIRCULARGRAPH Plot an interactive circular graph to illustrate connections in a network.
    %
    %% Syntax
    % circularGraph(X)
    % circularGraph(X,'PropertyName',propertyvalue,...)
    % h = circularGraph(...)
    %
    %% Description
    % A 'circular graph' is a visualization of a network of nodes and their
    % connections. The nodes are laid out along a circle, and the connections
    % are drawn within the circle. Click on a node to make the connections that
    % emanate from it more visible or less visible. Click on the 'Show All'
    % button to make all nodes and their connections visible. Click on the
    % 'Hide All' button to make all nodes and their connections less visible.
    %
    % Required input arguments.
    % X : A symmetric matrix of numeric or logical values.
    %
    % Optional properties.
    % Colormap : A N by 3 matrix of [r g b] triples, where N is the
    %            length(adjacenyMatrix).
    % Label    : A cell array of N strings.
    %%
    % Copyright 2016 The MathWorks, Inc.
    properties
        Node = node(0,0,1); % Array of nodes
        ColorMap;         % Colormap
        Label;            % Cell array of strings
        ShowButton;       % Turn all nodes on
        HideButton;       % Turn all nodes off
    end
    
    methods
        function this = circularGraph_noButtons_colorGradient(adjacencyMatrix,varargin)
            % Constructor
            p = inputParser;
            
            defaultColorMap = parula(length(adjacencyMatrix));
            defaultLabel = cell(length(adjacencyMatrix));
            for i = 1:length(defaultLabel)
                defaultLabel{i} = num2str(i);
            end
            
            addRequired(p,'adjacencyMatrix',@(x)(isnumeric(x) || islogical(x)));
            addParameter(p,'ColorMap',defaultColorMap,@(colormap)length(colormap) == 20);
            addParameter(p,'Label'   ,defaultLabel   ,@iscell);
            
            parse(p,adjacencyMatrix,varargin{:});
            this.ColorMap = p.Results.ColorMap;
            this.Label    = p.Results.Label;
            
%             this.ShowButton = uicontrol(...
%                 'Style','pushbutton',...
%                 'Position',[0 40 80 40],...
%                 'String','Show All',...
%                 'Callback',@circularGraph.showNodes,...
%                 'UserData',this);
            
%             this.HideButton = uicontrol(...
%                 'Style','pushbutton',...
%                 'Position',[0 0 80 40],...
%                 'String','Hide All',...
%                 'Callback',@circularGraph.hideNodes,...
%                 'UserData',this);
            
            fig = gcf;
            set(fig,...
                'UserData',this,...
                'CloseRequestFcn',@circularGraph_noButtons_nonLineScale.CloseRequestFcn);
            
     
            
            % Find non-zero values of s and their indices
            [row,col,v] = find(adjacencyMatrix);
            v = abs(v);
            % Calculate line widths based on values of s (stored in v).
            minLineWidth  = .3;
            lineWidthCoef = 3;%5;
            lineWidth = abs(v./5);%abs(v./max(v));
            if sum(lineWidth) == numel(lineWidth) % all lines are the same width.
                %lineWidth = repmat(3*minLineWidth,numel(lineWidth),1);
                lineWidth = lineWidthCoef*ones(length(lineWidth),1)/2 + minLineWidth;
            else % lines of variable width.
                lineWidth = lineWidthCoef*ones(length(lineWidth),1)/2 + minLineWidth;
            end
            
            rowSameAsCol = row == col;
            
            lineWidthNode = zeros(1, length(adjacencyMatrix)) + 0.001;
            lineWidthNode(row(rowSameAsCol)) = lineWidth(rowSameAsCol);
            
            v_range = 0:0.01:0.1;
            %v_range = [0 0.001 0.005 0.01 0.03 0.4 0.55 0.7 0.8 0.9 0.1];
            
            % Draw the nodes
            delete(this.Node);
            t = linspace(-pi,pi,length(adjacencyMatrix) + 1).'; % theta for each node
            extent = zeros(length(adjacencyMatrix),1);
            for i = 1:length(adjacencyMatrix)
                
                this_v = abs(adjacencyMatrix(row(i), col(i)));
                this_v
                if(this_v >= v_range(1) && this_v < v_range(2))
                    colorInd = 1;
                elseif(this_v >= v_range(2) && this_v < v_range(3))
                    colorInd = 2;
                elseif(this_v >= v_range(3) && this_v < v_range(4))
                    colorInd = 3;
                elseif(this_v >= v_range(4) && this_v < v_range(5))
                    colorInd = 4;    
                elseif(this_v >= v_range(5) && this_v < v_range(6))
                    colorInd = 5;
                elseif(this_v >= v_range(6) && this_v < v_range(7))
                    colorInd = 6;
                elseif(this_v >= v_range(7) && this_v < v_range(8))
                    colorInd = 7;
                elseif(this_v >= v_range(8) && this_v < v_range(9))
                    colorInd = 8;
                elseif(this_v >= v_range(9) && this_v < v_range(10))
                    colorInd = 9;
                elseif(this_v >= v_range(10) && this_v < v_range(11))
                    colorInd = 10;
                end
                this.Node(i) = node(cos(t(i)),sin(t(i)),lineWidthNode(i));
                if(lineWidthNode(i) == 0.001)
                    this.Node(i).Color = [1 1 1];
                elseif(adjacencyMatrix(i,i) == 0)
                    this.Node(i).Color = [1 1 1];
                elseif(adjacencyMatrix(i,i) > 0)
                    this.Node(i).Color = this.ColorMap(colorInd,:);
                elseif(adjacencyMatrix(i,i) < 0)
                    this.Node(i).Color = this.ColorMap(colorInd+10,:);
                elseif(isnan(adjacencyMatrix(i,i)))
                    
                end
                this.Node(i).Label = this.Label{i};
            end
            
            
            
            
            % Draw connections on the Poincare hyperbolic disk.
            %
            % Equation of the circles on the disk:
            % x^2 + y^2
            % + 2*(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1))*x
            % - 2*(u(1)-v(1))/(u(1)*v(2)-u(2)*v(1))*y + 1 = 0,
            % where u and v are points on the boundary.
            %
            % Standard form of equation of a circle
            % (x - x0)^2 + (y - y0)^2 = r^2
            %
            % Therefore we can identify
            % x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
            % y0 = (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
            % r^2 = x0^2 + y0^2 - 1
            
            for i = 1:length(v)
                this_v = abs(adjacencyMatrix(row(i), col(i)));
                if(this_v >= v_range(1) && this_v < v_range(2))
                    colorInd = 1;
                elseif(this_v >= v_range(2) && this_v < v_range(3))
                    colorInd = 2;
                elseif(this_v >= v_range(3) && this_v < v_range(4))
                    colorInd = 3;
                elseif(this_v >= v_range(4) && this_v < v_range(5))
                    colorInd = 4;    
                elseif(this_v >= v_range(5) && this_v < v_range(6))
                    colorInd = 5;
                elseif(this_v >= v_range(6) && this_v < v_range(7))
                    colorInd = 6;
                elseif(this_v >= v_range(7) && this_v < v_range(8))
                    colorInd = 7;
                elseif(this_v >= v_range(8) && this_v < v_range(9))
                    colorInd = 8;
                elseif(this_v >= v_range(9) && this_v < v_range(10))
                    colorInd = 9;
                elseif(this_v >= v_range(10) && this_v < v_range(11))
                    colorInd = 10;
                end
                
                if row(i) ~= col(i)
                    if abs(row(i) - col(i)) - length(adjacencyMatrix)/2 == 0
                        % points are diametric, so draw a straight line
                        u = [cos(t(row(i)));sin(t(row(i)))];
                        v = [cos(t(col(i)));sin(t(col(i)))];
                        
                        if(adjacencyMatrix(row(i), col(i)) > 0)
                            
                        elseif(adjacencyMatrix(row(i), col(i)) < 0)
                            colorInd = colorInd + 10;
                        end
                        
                        this.Node(row(i)).Connection(end+1) = line(...
                            [u(1);v(1)],...
                            [u(2);v(2)],...
                            'LineWidth', lineWidth(i),...
                            'Color', this.ColorMap(colorInd,:),...
                            'PickableParts','none');
                    else % points are not diametric, so draw an arc
                        u  = [cos(t(row(i)));sin(t(row(i)))];
                        v  = [cos(t(col(i)));sin(t(col(i)))];
                        x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
                        y0 =  (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
                        r  = sqrt(x0^2 + y0^2 - 1);
                        thetaLim(1) = atan2(u(2)-y0,u(1)-x0);
                        thetaLim(2) = atan2(v(2)-y0,v(1)-x0);
                        
                        if u(1) >= 0 && v(1) >= 0
                            % ensure the arc is within the unit disk
                            theta = [linspace(max(thetaLim),pi,50),...
                                linspace(-pi,min(thetaLim),50)].';
                        else
                            theta = linspace(thetaLim(1),thetaLim(2)).';
                        end
                        
                        i
                        this_v
                        if(adjacencyMatrix(row(i), col(i)) > 0)
                            
                        elseif(adjacencyMatrix(row(i), col(i)) < 0)
                            colorInd = colorInd + 10;
                        elseif(isnan(adjacencyMatrix(row(i), col(i))))
                            
                        end
                        
                        this.Node(row(i)).Connection(end+1) = line(...
                            r*cos(theta)+x0,...
                            r*sin(theta)+y0,...
                            'LineWidth', lineWidth(i),...
                            'Color', this.ColorMap(colorInd,:),...
                            'PickableParts','none');
                    end
                end
            end
            
            axis image;
            ax = gca;
            for i = 1:length(adjacencyMatrix)
                extent(i) = this.Node(i).Extent;
            end
            extent = max(extent(:));
            ax.XLim = ax.XLim + extent*[-1 1];
            fudgeFactor = 1.75; % Not sure why this is necessary. Eyeballed it.
            ax.YLim = ax.YLim + fudgeFactor*extent*[-1 1];
            ax.Visible = 'off';
            ax.SortMethod = 'depth';
            
            fig = gcf;
            fig.Color = [1 1 1];
        end
        
    end
    
    methods (Static = true)
        function showNodes(this,~)
            % Callback for 'Show All' button
            n = this.UserData.Node;
            for i = 1:length(n)
                n(i).Visible = true;
            end
        end
        
        function hideNodes(this,~)
            % Callback for 'Hide All' button
            n = this.UserData.Node;
            for i = 1:length(n)
                n(i).Visible = false;
            end
        end
        
        function CloseRequestFcn(this,~)
            % Callback for figure CloseRequestFcn
            c = this.UserData;
            for i = 1:length(c.Node)
                delete(c.Node(i));
            end
            delete(gcf);
        end
        
    end
    
end


%     % CIRCULARGRAPH Plot an interactive circular graph to illustrate connections in a network.
%     %
%     %% Syntax
%     % circularGraph(X)
%     % circularGraph(X,'PropertyName',propertyvalue,...)
%     % h = circularGraph(...)
%     %
%     %% Description
%     % A 'circular graph' is a visualization of a network of nodes and their
%     % connections. The nodes are laid out along a circle, and the connections
%     % are drawn within the circle. Click on a node to make the connections that
%     % emanate from it more visible or less visible. Click on the 'Show All'
%     % button to make all nodes and their connections visible. Click on the
%     % 'Hide All' button to make all nodes and their connections less visible.
%     %
%     % Required input arguments.
%     % X : A symmetric matrix of numeric or logical values.
%     %
%     % Optional properties.
%     % Colormap : A N by 3 matrix of [r g b] triples, where N is the
%     %            length(adjacenyMatrix).
%     % Label    : A cell array of N strings.
%     %%
%     % Copyright 2016 The MathWorks, Inc.
%     properties
%         Node = node(0,0); % Array of nodes
%         ColorMap;         % Colormap
%         Label;            % Cell array of strings
%         ShowButton;       % Turn all nodes on
%         HideButton;       % Turn all nodes off
%     end
%     
%     methods
%         function this = circularGraph_noButtons(adjacencyMatrix,varargin)
%             % Constructor
%             p = inputParser;
%             
%             defaultColorMap = parula(length(adjacencyMatrix));
%             defaultLabel = cell(length(adjacencyMatrix));
%             for i = 1:length(defaultLabel)
%                 defaultLabel{i} = num2str(i);
%             end
%             
%             addRequired(p,'adjacencyMatrix',@(x)(isnumeric(x) || islogical(x)));
%             addParameter(p,'ColorMap',defaultColorMap,@(colormap)length(colormap) == length(adjacencyMatrix));
%             addParameter(p,'Label'   ,defaultLabel   ,@iscell);
%             
%             parse(p,adjacencyMatrix,varargin{:});
%             this.ColorMap = p.Results.ColorMap;
%             this.Label    = p.Results.Label;
%             
% %             this.ShowButton = uicontrol(...
% %                 'Style','pushbutton',...
% %                 'Position',[0 40 80 40],...
% %                 'String','Show All',...
% %                 'Callback',@circularGraph.showNodes,...
% %                 'UserData',this);
%             
% %             this.HideButton = uicontrol(...
% %                 'Style','pushbutton',...
% %                 'Position',[0 0 80 40],...
% %                 'String','Hide All',...
% %                 'Callback',@circularGraph.hideNodes,...
% %                 'UserData',this);
%             
%             fig = gcf;
%             set(fig,...
%                 'UserData',this,...
%                 'CloseRequestFcn',@circularGraph_noButtons.CloseRequestFcn);
%             
%             % Draw the nodes
%             delete(this.Node);
%             t = linspace(-pi,pi,length(adjacencyMatrix) + 1).'; % theta for each node
%             extent = zeros(length(adjacencyMatrix),1);
%             for i = 1:length(adjacencyMatrix)
%                 this.Node(i) = node(cos(t(i)),sin(t(i)));
%                 this.Node(i).Color = this.ColorMap(i,:);
%                 this.Node(i).Label = this.Label{i};
%             end
%             
%             % Find non-zero values of s and their indices
%             [row,col,v] = find(adjacencyMatrix);
%             
%             % Calculate line widths based on values of s (stored in v).
%             minLineWidth  = .5;
%             lineWidthCoef = 5;
%             lineWidth = v./max(v);
%             if sum(lineWidth) == numel(lineWidth) % all lines are the same width.
%                 lineWidth = repmat(minLineWidth,numel(lineWidth),1);
%             else % lines of variable width.
%                 lineWidth = lineWidthCoef*lineWidth + minLineWidth;
%             end
%             
%             % Draw connections on the Poincare hyperbolic disk.
%             %
%             % Equation of the circles on the disk:
%             % x^2 + y^2
%             % + 2*(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1))*x
%             % - 2*(u(1)-v(1))/(u(1)*v(2)-u(2)*v(1))*y + 1 = 0,
%             % where u and v are points on the boundary.
%             %
%             % Standard form of equation of a circle
%             % (x - x0)^2 + (y - y0)^2 = r^2
%             %
%             % Therefore we can identify
%             % x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
%             % y0 = (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
%             % r^2 = x0^2 + y0^2 - 1
%             
%             for i = 1:length(v)
%                 if row(i) ~= col(i)
%                     if abs(row(i) - col(i)) - length(adjacencyMatrix)/2 == 0
%                         % points are diametric, so draw a straight line
%                         u = [cos(t(row(i)));sin(t(row(i)))];
%                         v = [cos(t(col(i)));sin(t(col(i)))];
%                         this.Node(row(i)).Connection(end+1) = line(...
%                             [u(1);v(1)],...
%                             [u(2);v(2)],...
%                             'LineWidth', lineWidth(i),...
%                             'Color', this.ColorMap(row(i),:),...
%                             'PickableParts','none');
%                     else % points are not diametric, so draw an arc
%                         u  = [cos(t(row(i)));sin(t(row(i)))];
%                         v  = [cos(t(col(i)));sin(t(col(i)))];
%                         x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
%                         y0 =  (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
%                         r  = sqrt(x0^2 + y0^2 - 1);
%                         thetaLim(1) = atan2(u(2)-y0,u(1)-x0);
%                         thetaLim(2) = atan2(v(2)-y0,v(1)-x0);
%                         
%                         if u(1) >= 0 && v(1) >= 0
%                             % ensure the arc is within the unit disk
%                             theta = [linspace(max(thetaLim),pi,50),...
%                                 linspace(-pi,min(thetaLim),50)].';
%                         else
%                             theta = linspace(thetaLim(1),thetaLim(2)).';
%                         end
%                         
%                         this.Node(row(i)).Connection(end+1) = line(...
%                             r*cos(theta)+x0,...
%                             r*sin(theta)+y0,...
%                             'LineWidth', lineWidth(i),...
%                             'Color', this.ColorMap(row(i),:),...
%                             'PickableParts','none');
%                     end
%                 end
%             end
%             
%             axis image;
%             ax = gca;
%             for i = 1:length(adjacencyMatrix)
%                 extent(i) = this.Node(i).Extent;
%             end
%             extent = max(extent(:));
%             ax.XLim = ax.XLim + extent*[-1 1];
%             fudgeFactor = 1.75; % Not sure why this is necessary. Eyeballed it.
%             ax.YLim = ax.YLim + fudgeFactor*extent*[-1 1];
%             ax.Visible = 'off';
%             ax.SortMethod = 'depth';
%             
%             fig = gcf;
%             fig.Color = [1 1 1];
%         end
%         
%     end
%     
%     methods (Static = true)
%         function showNodes(this,~)
%             % Callback for 'Show All' button
%             n = this.UserData.Node;
%             for i = 1:length(n)
%                 n(i).Visible = true;
%             end
%         end
%         
%         function hideNodes(this,~)
%             % Callback for 'Hide All' button
%             n = this.UserData.Node;
%             for i = 1:length(n)
%                 n(i).Visible = false;
%             end
%         end
%         
%         function CloseRequestFcn(this,~)
%             % Callback for figure CloseRequestFcn
%             c = this.UserData;
%             for i = 1:length(c.Node)
%                 delete(c.Node(i));
%             end
%             delete(gcf);
%         end
%         
%     end
%     
% end