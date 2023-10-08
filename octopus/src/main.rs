use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use regex::Regex;
use std::fmt;

pub struct SVEvent {
    chrom: String,
    pos: i32,
    id: String,
    ref: String,
    alt: String,
    qual: String,
    filter: String,
    info: HashMap<String, Option<String>>,
    format: String,
    sample: String,
}

impl SVEvent {
    fn new(chrom: String, pos: i32, id: String, ref: String, alt: String, qual: String, filter: String, info: String, format: String, sample: String) -> SVEvent {
        SVEvent {
            chrom,
            pos,
            id,
            ref,
            alt,
            qual,
            filter,
            info: SVEvent::parse_info(&info),
            format,
            sample,
        }
    }

    fn parse_info(info: &str) -> HashMap<String, Option<String>> {
        let mut info_dict = HashMap::new();
        for item in info.split(';') {
            let parts: Vec<&str> = item.split('=').collect();
            if parts.len() == 2 {
                info_dict.insert(parts[0].to_string(), Some(parts[1].to_string()));
            } else {
                info_dict.insert(parts[0].to_string(), None);
            }
        }
        info_dict
    }
} 





fn check_vcf_format(vcf_file_path: &str) -> io::Result<()> {
    let path = Path::new(vcf_file_path);
    let file = File::open(&path)?;
    let reader = io::BufReader::new(file);

    let mut has_header = false;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            has_header = true;
            continue;
        }

        if line.contains(' ') {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "ERROR: Invalid VCF format. Non-header lines should not contain spaces."));
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, format!("ERROR: Invalid VCF format. Expected at least 10 fields, but got {}", fields.len())));
        }

        if fields[1].parse::<i32>().is_err() {
            return Err(io::Error::new(io::ErrorKind::InvalidData, format!("ERROR: Invalid VCF format. Position (field 2) should be a number, but got {}", fields[1])));
        }

        if fields[5] != "." && fields[5].parse::<f32>().is_err() {
            return Err(io::Error::new(io::ErrorKind::InvalidData, format!("ERROR: Invalid VCF format. Quality score (field 6) should be a number or '.', but got {}", fields[5])));
        }
    }

    if !has_header {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "ERROR: Invalid VCF format. The file should contain at least one header line starting with '#'"));
    }

    Ok(())
}




struct Event {
    chrom: String,
    alt: String,
}

impl Event {
    fn is_bnd(&self) -> bool {
        // Implement this function based on your logic
        true
    }
}

fn is_same_chr_bnd(event: &Event) -> bool {
    if event.is_bnd() {
        let re = Regex::new(r"[\[\]:]").unwrap();
        let split_result: Vec<&str> = re.split(&event.alt).collect();
        if split_result.len() != 4 {
            println!("Unexpected ALT format, it should be something like N]chr10:69650962]: {:?}", split_result);
        } else {
            let chrom_alt = split_result[1];
            return event.chrom == chrom_alt;
        }
    }
    false  // For non-BND, we won't categorize them as same_chr_bnd or diff_chr_bnd events
}

fn get_bnd_pattern(alt: &str) -> Option<&str> {
    let first_char = alt.chars().nth(0).unwrap();
    let second_char = alt.chars().nth(1).unwrap();
    let last_char = alt.chars().last().unwrap();
    let second_last_char = alt.chars().nth(alt.len() - 2).unwrap();

    if "ATCGN".contains(first_char) && second_char == '[' {
        Some("t[p[")
    } else if "ATCGN".contains(first_char) && second_char == ']' {
        Some("t]p]")
    } else if "ATCGN".contains(last_char) && second_last_char == '[' {
        Some("[p[t")
    } else if "ATCGN".contains(last_char) && second_last_char == ']' {
        Some("]p]t")
    } else {
        None
    }
}


fn get_alt_chrom_pos(alt: &str) -> (Option<String>, Option<i32>) {
    let re = Regex::new(r"[\[\]:]").unwrap();
    let split_result: Vec<&str> = re.split(alt).collect();
    if split_result.len() != 4 {
        println!("Unexpected ALT format, it should be something like N]chr10:69650962]: {:?}", split_result);
        return (None, None);
    } else {
        let chrom_alt = split_result[1].to_string();
        let pos_alt: i32 = split_result[2].parse().unwrap();
        return (Some(chrom_alt), Some(pos_alt));
    }
}

fn find_mate_bnd_and_no_mate_events(events: Vec<Event>, pos_tolerance: i32) -> (Vec<(Event, Event)>, Vec<Event>) {
    let mut event_dict: HashMap<(String, i32, String, i32), Event> = HashMap::new();
    let mut mate_bnd_pairs: Vec<(Event, Event)> = Vec::new();
    let mut no_mate_events: Vec<Event> = Vec::new();

    for event in events {
        let (chrom_alt, pos_alt) = get_alt_chrom_pos(&event.alt);
        let key = (event.chrom.clone(), event.pos, chrom_alt.unwrap(), pos_alt.unwrap());

        let mut mate_found = false;
        for i in -pos_tolerance..=pos_tolerance {
            for j in -pos_tolerance..=pos_tolerance {
                let reverse_key = (chrom_alt.clone().unwrap(), pos_alt.unwrap() + i, event.chrom.clone(), event.pos + j);
                if event_dict.contains_key(&reverse_key) {
                    let mate_event = event_dict.remove(&reverse_key).unwrap();
                    mate_bnd_pairs.push((mate_event, event.clone()));
                    mate_found = true;
                    break;
                }
            }
            if mate_found {
                break;
            }
        }

        if !mate_found {
            event_dict.insert(key, event);
        }
    }

    no_mate_events = event_dict.values().cloned().collect();
    return (mate_bnd_pairs, no_mate_events);
}












fn get_bnd_pattern(alt: &str) -> Option<&str> {
    // Get the pattern of BND: t[p[, t]p], ]p]t, [p[t
    match alt.chars().nth(1) {
        Some('[') => Some("t[p["),
        Some(']') => Some("t]p]"),
        _ => None,
    }
}

fn is_mate_pair_reciprocal_translocation(event1: &Event, event2: &Event) -> bool {
    let qualified_pairings = vec![
        vec!["]p]t", "]p]t"],
        vec!["t[p[", "t[p["],
        vec!["t]p]", "[p[t"],
    ];
    let event1_pattern = get_bnd_pattern(&event1.alt);
    let event2_pattern = get_bnd_pattern(&event2.alt);

    qualified_pairings.into_iter().any(|pairing| {
        vec![event1_pattern, event2_pattern].sort() == pairing.sort()
    })
}

fn is_same_bnd_event(event1: &Event, event2: &Event) -> bool {
    let qualified_pairings = vec![
        vec!["]p]t", "t[p["]
    ];
    let event1_pattern = get_bnd_pattern(&event1.alt);
    let event2_pattern = get_bnd_pattern(&event2.alt);

    qualified_pairings.into_iter().any(|pairing| {
        vec![event1_pattern, event2_pattern].sort() == pairing.sort()
    })
}

fn compare_chromosomes(event1: &Event, event2: &Event) -> bool {
    let original_order = vec![&event1.chrom, &event2.chrom];
    let mut sorted_order = original_order.clone();
    sorted_order.sort();

    original_order == sorted_order
}

fn is_independent_bnd_event(event1: &Event, event2: &Event) -> bool {
    let qualified_pairings = vec![
        vec!["t]p]", "]p]t"],
        vec!["t]p]", "t]p]"],
        vec!["[p[t", "[p[t"],
        vec!["t[p[", "t]p]"],
        vec!["t[p[", "[p[t"],
        vec!["]p]t", "[p[t"]
    ];
    let event1_pattern = get_bnd_pattern(&event1.alt);
    let event2_pattern = get_bnd_pattern(&event2.alt);

    qualified_pairings.into_iter().any(|pairing| {
        vec![event1_pattern, event2_pattern].sort() == pairing.sort()
    })
}






fn parse_vcf(vcf_file_path: &str) -> io::Result<(Vec<String>, Vec<SVEvent>, Vec<SVEvent>, Vec<SVEvent>)> {
    let mut same_chr_bnd_events = Vec::new();
    let mut diff_chr_bnd_events = Vec::new();
    let mut non_bnd_events = Vec::new();
    let mut headers = Vec::new();

    let file = File::open(&Path::new(vcf_file_path))?;
    for line in io::BufReader::new(file).lines() {
        let line = line?;
        if line.starts_with('#') {
            headers.push(line);
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        let mut info = HashMap::new();
        for item in fields[7].split(';') {
            let parts: Vec<&str> = item.split('=').collect();
            if parts.len() == 2 {
                info.insert(parts[0].to_string(), parts[1].to_string());
            } else {
                info.insert(parts[0].to_string(), "".to_string());
            }
        }

        let event = SVEvent::new(
            fields[0].to_string(),
            fields[1].parse().unwrap(),
            fields[2].to_string(),
            fields[3].to_string(),
            fields[4].to_string(),
            fields[5].to_string(),
            fields[6].to_string(),
            info,
            fields[8].to_string(),
            fields[9].to_string(),
        );

        if event.is_bnd() {
            if is_same_chr_bnd(&event) {
                same_chr_bnd_events.push(event);
            } else {
                diff_chr_bnd_events.push(event);
            }
        } else {
            non_bnd_events.push(event);
        }
    }

    Ok((headers, same_chr_bnd_events, diff_chr_bnd_events, non_bnd_events))
}

fn is_same_chr_bnd(event: &SVEvent) -> bool {
    // Implement this function
    unimplemented!()
}

fn get_alt_chrom_pos(alt: &str) -> (Option<String>, Option<i32>) {
    // Implement this function
    unimplemented!()
}

fn find_special_no_mate_diff_bnd_pair_and_other_single_TRA(events: Vec<SVEvent>, pos_tolerance: i32) -> (Vec<(SVEvent, SVEvent)>, Vec<SVEvent>) {
    // Implement this function
    unimplemented!()
}

fn is_special_no_mate_diff_bnd_pair_reciprocal_translocation(event1: &SVEvent, event2: &SVEvent) -> bool {
    // Implement this function
    unimplemented!()
}

fn is_independent_special_bnd_event(event1: &SVEvent, event2: &SVEvent) -> bool {
    // Implement this function
    unimplemented!()
}




fn main() {
    println!("Hello, world!");
}
