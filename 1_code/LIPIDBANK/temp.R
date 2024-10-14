
all_word2 <-
  all_word[-grep("bacteria|human|infant|patient|dog|rat|vertebrates|rabbit|snake|chiken|fungal|cattle|kangaroos|opossum|pig|koala|mammalian|neonatal|python|fungus|fungu|pelicans",
                 stringr::str_to_lower(all_word))] %>%
  stringr::str_to_lower()

all_word2 <-
  all_word2[-grep("owls|bean|newborn|baboon|kite|serum|urine|rana|bullfrog|bufo|varanus|amyda|fish",
                  all_word2)]

all_word2 <-
  all_word2[-grep("moth|peach|drasche|ginseng|sertifer|cyperus|asteraceae|pine|bollworm|sheep",
                  all_word2)]

all_word2 <-
  all_word2[-grep("carrot|wheat|milk|lemon|cabbage|sulfolobus|shark|methanosphaera",
                  all_word2)]

all_word2 <-
  all_word2[-grep("barin|lung|spleen|chicken|retina|alga|lobster|pyrococcus|caldariella|mucor|mouse",
                  all_word2)]

all_word2 <-
  all_word2[-grep("paca|moschatus|monkey|marronnier|sponge|plant|bacterium|brevundimonas|codium",
                  all_word2)]

all_word2 <-
  all_word2[-grep(" sp",
                  all_word2)]

all_word2 <-
  all_word2[-grep("boar|cortex|bovine|rhamnus|cinchona|cotton|female|trillium|dioscorea|seed",
                  all_word2)]

all_word2 <-
  all_word2[-grep("digitonin|diginatin|root|sarsaparilla|scilla|water|porcine|pregnant|proteus|chick|pseudomonas",
                  all_word2)]

all_word2 <-
  all_word2[-grep("rhizobium|plasma|corn|soil|ruscus|petals|rodent|ricinus|rice|rhodovulum",
                  all_word2)]

all_word2 <-
  all_word2[-grep("rhodospirillum|rhodocyclus|mytilus|prasinophyceae|sansevieria|salmonella",
                  all_word2)]

all_word2 <-
  all_word2[-grep("schizosaccharomyces|schizophylum|scallop|hydnocarpus|quinqueradiata|shigella|worm",
                  all_word2)]

all_word2 <-
  all_word2[-grep("primates|palm|acinetobacter|acnistus|actinobacillus|aeromonas|aeropyrum|alligator|amaroucium|tunicate",
                  all_word2)]

all_word2 <-
  all_word2[-grep("amphiuma|animal|bird|blood|bordetella|brain|caiman|campylobacter|capsicum|feces",
                  all_word2)]

all_word2 <-
  all_word2[-grep("citrullus|cocoa|codiaceae|yersinia|yeast|whale|vibrio|liver|vibrio|scammony|tissue|fruit",
                  all_word2)]

all_word2 <-
  all_word2[-grep("ranunculaceae|xenopus|xanthomonas|malvacearum|xylininum|bee|grease|wool|solanaceae|urinary|wax",
                  all_word2)]

all_word2 <-
  all_word2[-grep("latirostris|sultana|visceral|violet|vinegar|ustilago|liliaceae|uredovora|viscera|trichosanthes",
                  all_word2)]

all_word2 <-
  all_word2[-grep("ant|tragopogon|ascidiacea|marine|achlya|acetobacter|grape|asteroidea|aspergillus|asclepias|apple",
                  all_word2)]

all_word2 <-
  all_word2[-grep("umbellatum|tanghinin|clam|berries|tanghinin",
                  all_word2)]

all_word2 <-
  all_word2[-grep("umbellatum|tanghinin|clam|berries|tanghinin|tabacco|syringae|creeper",
                  all_word2)]

all_word2 <-
  all_word2[-grep("synthetic product",
                  all_word2)]

all_word2 <-
  all_word2[-grep("prodaction mechanism",
                  all_word2)]

all_word2 <-
  all_word2[-grep("ox ",
                  all_word2)]

all_word2 <-
  all_word2[-grep("sea|sardine|pterosperma|providencia|porphyromonas|porifera|trilinolenin|butter|gonyaulax",
                  all_word2)]

all_word2 <-
  all_word2[-grep("toad|tilapia|thunnus|thevetia|thermococcus|thalictrum|telesto|synthetic|sunflower|sturgeon",
                  all_word2)]

all_word2 <-
  all_word2[-grep("streptococcus|strain|bacilli|coral|pernyi|tubercle|adrenal|alestes|leave|gestation|amoebae|anacystis",
                  all_word2)]

all_word2 <-
  all_word2[-grep("testis|apium|arapaima|skin|tree|banana|leaves|bacteroides|azospirillum|flower|rhodococcus",
                  all_word2)]

all_word2 <-
  all_word2[-grep("rhodobacter|rhizomes|rambutan|racemosa|pyramimonas|potamon|peucedanum|penicillium|pectinatus",
                  all_word2)]

all_word2 <-
  all_word2[-grep("pectinatus|pasteurianus|crab|oyster|oxkidney|potato|oncoryhnchus|olive|nocardia|neurospora",
                  all_word2)]

all_word2 <-
  all_word2[-grep("neriifolia|neisseria",
                  all_word2)]


all_word2 <-
  all_word2[!is.na(all_word2)]
length(all_word2)
sort(all_word2)

openxlsx::write.xlsx(data.frame(all_word = sort(all_word2)), file = "all_word.xlsx", asTable = TRUE)
