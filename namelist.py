from datetime import datetime, timedelta

col1_width = 20
col2_width = 25

# Format each line with fixed column widths
def format_line(col1, col2):
    return f"{col1:<{col1_width}}{col2:<{col2_width}}"


def generate_namelist(expno, start_str, start_time, end_str, restart, endmin, timeinterval):
    
    n = int(endmin / timeinterval)

    output = []
    output.append(format_line("1", "input namelist"))
    output.append(format_line(restart, "restart"))
    output.append(format_line(expno, "exp version"))
    output.append(format_line(start_str, "starting date"))
    output.append(format_line(end_str, "end date"))
    
    
    start = datetime(*start_time)
    
    intervals = [timedelta(minutes=timeinterval)]*n
    
    for i in range(n+1):
        output.append(format_line((start + sum(intervals[:i], timedelta())).strftime('%Y-%m-%d:%H:%M:%S'), "posting date"))  

    output.append(format_line("stop", "end of post date string"))

    return "\n".join(output)

expno = 74

output = generate_namelist(expno, '1990-01-01:00:00:00', [1990, 1, 1, 0, 0, 0], '1990-01-10:00:00:00', 0, 10*(24*60), 120)

with open("input_panarctic_muphi", "w") as file:
    file.write(output)